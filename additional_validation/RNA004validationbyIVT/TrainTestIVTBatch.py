import os
import csv
import hashlib
import numpy as np
import pandas as pd
from pyarrow import parquet as pq
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import ModelCheckpoint
import matplotlib.pyplot as plt
import nnmodels.CNNWavenet as cnnwavenet

import nnmodels.Transformer1D as transformer

INPUT_PARQUET = "/mnt/share/ueda/tycoon_revise/rna004_training/rna004_ivt_training.parquet"
# OUTPUT_DIR = "/mnt/share/ueda/tycoon_revise/rna004_cnnwavenet_test"
OUTPUT_DIR = "/mnt/share/ueda/tycoon_revise/rna004_transformer_test"

N_TRAIN_PER_CLASS = 10000
N_VALIDATION_PER_CLASS = 2000
EPOCHS = 100
BATCH_SIZE = 256
LEARNING_RATE = 0.0008
SEED = 100

# Force TensorFlow to use CPU only.
# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

def set_seed(seed):
    np.random.seed(seed)
    tf.random.set_seed(seed)


def stable_hash(text):
    return int(hashlib.md5(text.encode("utf-8")).hexdigest()[:8], 16)


def load_dataframe(path):
    table = pq.read_table(path, columns=["read_id", "trna", "trimsignal"])
    return table.to_pandas()


def assign_read_split(df):
    # Split at read level to avoid leakage between tRNA segments
    # originating from the same tandem IVT molecule.
    split = {}
    for read_id in sorted(df["read_id"].unique()):
        split[read_id] = "validation" if stable_hash(read_id) % 6 == 0 else "train"
    df = df.copy()
    df["dataset_split"] = df["read_id"].map(split)
    return df


def sample_per_class(
    df,
    trna_names,
):
    train_parts = []
    val_parts = []

    for trna in trna_names:

        sub = df[
            df["trna"] == trna
        ].copy()

        print(
            trna,
            "available:",
            len(sub),
        )

        required = (
            N_TRAIN_PER_CLASS
            + N_VALIDATION_PER_CLASS
        )

        if len(sub) < required:
            raise RuntimeError(
                "Not enough reads for %s: %d"
                % (
                    trna,
                    len(sub),
                )
            )

        # Use only the first 12,000 signals.
        sub = sub.iloc[
            :required
        ].copy()

        # First 10,000 for training.
        train = sub.iloc[
            :N_TRAIN_PER_CLASS
        ].copy()

        # Remaining 2,000 for validation.
        val = sub.iloc[
            N_TRAIN_PER_CLASS:
            N_TRAIN_PER_CLASS
            + N_VALIDATION_PER_CLASS
        ].copy()

        train_parts.append(
            train
        )

        val_parts.append(
            val
        )

        print(
            "  train:",
            len(train),
            "validation:",
            len(val),
        )

    train_df = pd.concat(
        train_parts,
        ignore_index=True,
    )

    val_df = pd.concat(
        val_parts,
        ignore_index=True,
    )

    return (
        train_df,
        val_df,
    )

def dataframe_to_arrays(df, trna_to_index):
    x = np.stack([np.asarray(v, dtype=np.float32) for v in df["trimsignal"]], axis=0)
    y_index = np.asarray([trna_to_index[v] for v in df["trna"]], dtype=np.int32)
    x = np.expand_dims(x, axis=-1)
    y = keras.utils.to_categorical(y_index, num_classes=len(trna_to_index))
    return x, y, y_index


def load_data(path):
    df = load_dataframe(path)
    trna_names = sorted(df["trna"].unique())
    if len(trna_names) != 16:
        raise RuntimeError("Expected 16 tRNA classes, found %d: %s" % (len(trna_names), trna_names))

    df = assign_read_split(df)
    train_df, val_df = sample_per_class(df, trna_names)
    trna_to_index = {trna: i for i, trna in enumerate(trna_names)}

    x_train, y_train, _ = dataframe_to_arrays(train_df, trna_to_index)
    x_val, y_val, y_val_index = dataframe_to_arrays(val_df, trna_to_index)

    print("Training:", x_train.shape, y_train.shape)
    print("Validation:", x_val.shape, y_val.shape)
    return x_train, y_train, x_val, y_val, y_val_index, trna_names


# def build_model(signal_length, num_classes):
#     shape = (None, signal_length, 1)
#     return cnnwavenet.build_network(shape=shape, num_classes=num_classes, do_r=0.2)

def build_model(
    signal_length,
    num_classes,
):
    shape = (
        None,
        signal_length,
        1,
    )

    return transformer.build_network(
        shape=shape,
        num_classes=num_classes,
        do_r=0.2,
    )

def train_model(x_train, y_train, x_val, y_val):
    set_seed(SEED)
    model = build_model(x_train.shape[1], y_train.shape[1])
    model.summary()

    optimizer = keras.optimizers.Adam(
        learning_rate=LEARNING_RATE,
        beta_1=0.9,
        beta_2=0.999,
    )

    model.compile(
        loss="categorical_crossentropy",
        optimizer=optimizer,
        metrics=["accuracy"],
    )

    weight_path = os.path.join(OUTPUT_DIR, "best_weights.h5")
    checkpoint = ModelCheckpoint(
        filepath=weight_path,
        monitor="val_accuracy",
        verbose=1,
        save_best_only=True,
        save_weights_only=True,
        mode="max",
    )

    history = model.fit(
        x_train,
        y_train,
        epochs=EPOCHS,
        batch_size=BATCH_SIZE,
        shuffle=True,
        verbose=1,
        validation_data=(x_val, y_val),
        callbacks=[checkpoint],
    )

    history_df = pd.DataFrame(history.history)
    history_df.to_csv(os.path.join(OUTPUT_DIR, "history.csv"), index=False)
    model.load_weights(weight_path)
    return model, history_df


def make_confusion_matrix(y_true, y_pred, num_classes):
    matrix = np.zeros((num_classes, num_classes), dtype=np.int64)
    for t, p in zip(y_true, y_pred):
        matrix[t, p] += 1
    return matrix


def save_confusion_matrix(matrix, trna_names):
    raw_df = pd.DataFrame(matrix, index=trna_names, columns=trna_names)
    raw_df.to_csv(os.path.join(OUTPUT_DIR, "confusion_matrix_counts.tsv"), sep="\t")

    denom = matrix.sum(axis=1, keepdims=True).astype(float)
    denom[denom == 0] = 1.0
    norm = matrix / denom
    pd.DataFrame(norm, index=trna_names, columns=trna_names).to_csv(
        os.path.join(OUTPUT_DIR, "confusion_matrix_normalized.tsv"), sep="\t"
    )

    plt.figure(figsize=(10, 9))
    plt.imshow(norm, aspect="auto", vmin=0.0, vmax=1.0)
    plt.colorbar(label="Fraction")
    plt.xticks(np.arange(len(trna_names)), trna_names, rotation=90)
    plt.yticks(np.arange(len(trna_names)), trna_names)
    plt.xlabel("Predicted tRNA")
    plt.ylabel("True tRNA")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "confusion_matrix.png"), dpi=300)
    plt.savefig(os.path.join(OUTPUT_DIR, "confusion_matrix.pdf"))
    plt.close()


def plot_learning_curves(history_df):
    epochs = np.arange(len(history_df)) + 1

    plt.figure(figsize=(7, 5))
    plt.plot(epochs, history_df["accuracy"], label="Training", linewidth=1.2)
    plt.plot(epochs, history_df["val_accuracy"], label="Validation", linewidth=1.2, linestyle="--")
    plt.xlabel("Epoch")
    plt.ylabel("Accuracy")
    plt.xlim(1, len(history_df))
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "learning_curve_accuracy.png"), dpi=300)
    plt.savefig(os.path.join(OUTPUT_DIR, "learning_curve_accuracy.pdf"))
    plt.close()

    plt.figure(figsize=(7, 5))
    plt.plot(epochs, history_df["loss"], label="Training", linewidth=1.2)
    plt.plot(epochs, history_df["val_loss"], label="Validation", linewidth=1.2, linestyle="--")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.xlim(1, len(history_df))
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "learning_curve_loss.png"), dpi=300)
    plt.savefig(os.path.join(OUTPUT_DIR, "learning_curve_loss.pdf"))
    plt.close()


def save_class_index(trna_names):
    with open(os.path.join(OUTPUT_DIR, "tRNAindex.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["index", "trna"])
        for i, trna in enumerate(trna_names):
            writer.writerow([i, trna])


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    x_train, y_train, x_val, y_val, y_val_index, trna_names = load_data(INPUT_PARQUET)
    save_class_index(trna_names)

    model, history_df = train_model(x_train, y_train, x_val, y_val)
    prob = model.predict(x_val, batch_size=BATCH_SIZE, verbose=1)
    y_pred = np.argmax(prob, axis=1)

    matrix = make_confusion_matrix(y_val_index, y_pred, len(trna_names))
    save_confusion_matrix(matrix, trna_names)
    plot_learning_curves(history_df)

    accuracy = float(np.mean(y_pred == y_val_index))
    best_epoch = int(history_df["val_accuracy"].idxmax() + 1)
    best_val_accuracy = float(history_df["val_accuracy"].max())

    pd.DataFrame([{
        "validation_accuracy": accuracy,
        "best_epoch": best_epoch,
        "best_val_accuracy": best_val_accuracy,
        "n_train_per_class": N_TRAIN_PER_CLASS,
        "n_validation_per_class": N_VALIDATION_PER_CLASS,
        "epochs": EPOCHS,
        "batch_size": BATCH_SIZE,
        "learning_rate": LEARNING_RATE,
    }]).to_csv(os.path.join(OUTPUT_DIR, "summary.tsv"), sep="\t", index=False)

    print("Validation accuracy:", accuracy)
    print("Best epoch:", best_epoch)
    print("Best validation accuracy:", best_val_accuracy)


if __name__ == "__main__":
    main()
