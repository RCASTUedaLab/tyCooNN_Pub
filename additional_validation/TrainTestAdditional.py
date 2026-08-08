import os

# Force TensorFlow to use CPU only.
# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import csv
import numpy as np
import pandas as pd

from pyarrow import parquet as pq

import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.callbacks import ModelCheckpoint

import matplotlib.pyplot as plt

import nnmodels.CNNOnly as cnnonly
import nnmodels.CNNWavenet as cnnwavenet


# ============================================================
# Configuration
# ============================================================

INPUT_PARQUET = (
    "/mnt/share/ueda/tycoon_revise/"
    "rna004_training/"
    "rna004_ivt_training.parquet"
)

OUTPUT_DIR = (
    "/mnt/share/ueda/tycoon_revise/"
    "rna004_architecture_comparison"
)

MODEL_NAMES = [
    "CNNOnly",
    "CNNWavenet",
]

N_TRAIN_PER_CLASS = 10000
N_VALIDATION_PER_CLASS = 2000

EPOCHS = 100
BATCH_SIZE = 256
LEARNING_RATE = 0.0008

SEED = 100


# ============================================================
# Random seed
# ============================================================

def set_seed(seed):
    """
    Set NumPy and TensorFlow random seeds.
    """

    np.random.seed(
        seed
    )

    tf.random.set_seed(
        seed
    )


# ============================================================
# Load parquet
# ============================================================

def load_dataframe(path):
    """
    Load read ID, tRNA label, and fixed-length signal
    from the training parquet file.
    """

    table = pq.read_table(
        path,
        columns=[
            "read_id",
            "trna",
            "trimsignal",
        ],
    )

    return table.to_pandas()


# ============================================================
# Exact train / validation split
# ============================================================

def split_per_class(
    df,
    trna_names,
):
    """
    Use exactly 12,000 signals per tRNA class.

    For each class:
        10,000 signals -> training
         2,000 signals -> validation

    The split is deterministic because a fixed random seed
    is used.
    """

    train_parts = []
    validation_parts = []

    rng = np.random.RandomState(
        SEED
    )

    print()
    print(
        "=============================================="
    )

    print(
        "Creating exact train / validation split"
    )

    print(
        "=============================================="
    )

    print()

    for trna in trna_names:

        sub = df[
            df[
                "trna"
            ] == trna
        ].copy()

        n_available = len(
            sub
        )

        n_required = (
            N_TRAIN_PER_CLASS
            + N_VALIDATION_PER_CLASS
        )

        if (
            n_available
            < n_required
        ):

            raise RuntimeError(
                "Not enough signals for %s: "
                "available=%d required=%d"
                % (
                    trna,
                    n_available,
                    n_required,
                )
            )

        indices = sub.index.values.copy()

        rng.shuffle(
            indices
        )

        # Use exactly 12,000 signals.
        indices = indices[
            :n_required
        ]

        train_indices = indices[
            :N_TRAIN_PER_CLASS
        ]

        validation_indices = indices[
            N_TRAIN_PER_CLASS:
            N_TRAIN_PER_CLASS
            + N_VALIDATION_PER_CLASS
        ]

        train = sub.loc[
            train_indices
        ].copy()

        validation = sub.loc[
            validation_indices
        ].copy()

        train[
            "dataset_split"
        ] = "train"

        validation[
            "dataset_split"
        ] = "validation"

        train_parts.append(
            train
        )

        validation_parts.append(
            validation
        )

        print(
            "%-8s train=%5d validation=%4d total=%5d"
            % (
                trna,
                len(
                    train
                ),
                len(
                    validation
                ),
                len(
                    train
                )
                + len(
                    validation
                ),
            )
        )

    train_df = pd.concat(
        train_parts,
        ignore_index=True,
    )

    validation_df = pd.concat(
        validation_parts,
        ignore_index=True,
    )

    print()

    print(
        "Total training signals:",
        len(
            train_df
        ),
    )

    print(
        "Total validation signals:",
        len(
            validation_df
        ),
    )

    return (
        train_df,
        validation_df,
    )


# ============================================================
# Save exact dataset split
# ============================================================

def save_dataset_split(
    train_df,
    validation_df,
):
    """
    Save read IDs and class labels used for the exact split.
    This allows the same split to be reused later.
    """

    split_df = pd.concat(
        [
            train_df[
                [
                    "read_id",
                    "trna",
                    "dataset_split",
                ]
            ],
            validation_df[
                [
                    "read_id",
                    "trna",
                    "dataset_split",
                ]
            ],
        ],
        ignore_index=True,
    )

    output_path = os.path.join(
        OUTPUT_DIR,
        "dataset_split.tsv",
    )

    split_df.to_csv(
        output_path,
        sep="\t",
        index=False,
    )

    print(
        "Saved dataset split:",
        output_path,
    )


# ============================================================
# Convert dataframe to NumPy arrays
# ============================================================

def dataframe_to_arrays(
    df,
    trna_to_index,
):
    """
    Convert fixed-length signals and tRNA labels
    to arrays suitable for Keras.
    """

    x = np.stack(
        [
            np.asarray(
                signal,
                dtype=np.float32,
            )
            for signal
            in df[
                "trimsignal"
            ]
        ],
        axis=0,
    )

    y_index = np.asarray(
        [
            trna_to_index[
                trna
            ]
            for trna
            in df[
                "trna"
            ]
        ],
        dtype=np.int32,
    )

    # Add channel dimension.
    x = np.expand_dims(
        x,
        axis=-1,
    )

    y = keras.utils.to_categorical(
        y_index,
        num_classes=len(
            trna_to_index
        ),
    )

    return (
        x,
        y,
        y_index,
    )


# ============================================================
# Load and prepare dataset
# ============================================================

def load_data(
    path,
):
    """
    Load parquet and create an exact balanced split.
    """

    print()
    print(
        "Loading parquet:"
    )

    print(
        path
    )

    df = load_dataframe(
        path
    )

    trna_names = sorted(
        df[
            "trna"
        ].unique()
    )

    if (
        len(
            trna_names
        )
        != 16
    ):

        raise RuntimeError(
            "Expected 16 tRNA classes, found %d: %s"
            % (
                len(
                    trna_names
                ),
                trna_names,
            )
        )

    print()
    print(
        "tRNA classes:"
    )

    for i, trna in enumerate(
        trna_names
    ):

        print(
            "%2d %s"
            % (
                i,
                trna,
            )
        )

    (
        train_df,
        validation_df,
    ) = split_per_class(
        df,
        trna_names,
    )

    save_dataset_split(
        train_df,
        validation_df,
    )

    trna_to_index = {
        trna: i
        for i, trna
        in enumerate(
            trna_names
        )
    }

    (
        x_train,
        y_train,
        y_train_index,
    ) = dataframe_to_arrays(
        train_df,
        trna_to_index,
    )

    (
        x_validation,
        y_validation,
        y_validation_index,
    ) = dataframe_to_arrays(
        validation_df,
        trna_to_index,
    )

    print()
    print(
        "Training shape:",
        x_train.shape,
        y_train.shape,
    )

    print(
        "Validation shape:",
        x_validation.shape,
        y_validation.shape,
    )

    print(
        "Signal length:",
        x_train.shape[
            1
        ],
    )

    return (
        x_train,
        y_train,
        x_validation,
        y_validation,
        y_validation_index,
        trna_names,
    )


# ============================================================
# Save class index
# ============================================================

def save_class_index(
    trna_names,
):
    """
    Save class index used by all models.
    """

    output_path = os.path.join(
        OUTPUT_DIR,
        "tRNAindex.csv",
    )

    with open(
        output_path,
        "w",
        newline="",
    ) as handle:

        writer = csv.writer(
            handle
        )

        writer.writerow(
            [
                "index",
                "trna",
            ]
        )

        for i, trna in enumerate(
            trna_names
        ):

            writer.writerow(
                [
                    i,
                    trna,
                ]
            )


# ============================================================
# Build model
# ============================================================

def build_model(
    model_name,
    signal_length,
    num_classes,
):
    """
    Build CNN-only or CNN-WaveNet model.
    """

    shape = (
        None,
        signal_length,
        1,
    )

    if (
        model_name
        == "CNNOnly"
    ):

        model = cnnonly.build_network(
            shape=shape,
            num_classes=num_classes,
            do_r=0.2,
        )

    elif (
        model_name
        == "CNNWavenet"
    ):

        model = cnnwavenet.build_network(
            shape=shape,
            num_classes=num_classes,
            do_r=0.2,
        )

    else:

        raise ValueError(
            "Unknown model: %s"
            % model_name
        )

    return model


# ============================================================
# Train model
# ============================================================

def train_model(
    model_name,
    x_train,
    y_train,
    x_validation,
    y_validation,
    output_dir,
):
    """
    Train one architecture for 100 epochs.

    The checkpoint with the highest validation accuracy
    is retained and reloaded after training.
    """

    print()
    print(
        "=================================================="
    )

    print(
        "Training model:",
        model_name
    )

    print(
        "=================================================="
    )

    set_seed(
        SEED
    )

    model = build_model(
        model_name=model_name,
        signal_length=x_train.shape[
            1
        ],
        num_classes=y_train.shape[
            1
        ],
    )

    model.summary()

    print()
    print(
        "Model output shape:",
        model.output_shape,
    )

    if (
        model.output_shape[
            -1
        ]
        != 16
    ):

        raise RuntimeError(
            "Model output layer does not contain 16 classes: %s"
            % (
                model.output_shape,
            )
        )

    optimizer = keras.optimizers.Adam(
        learning_rate=LEARNING_RATE,
        beta_1=0.9,
        beta_2=0.999,
    )

    model.compile(
        loss="categorical_crossentropy",
        optimizer=optimizer,
        metrics=[
            "accuracy"
        ],
    )

    weight_path = os.path.join(
        output_dir,
        "best_weights.h5",
    )

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
        validation_data=(
            x_validation,
            y_validation,
        ),
        callbacks=[
            checkpoint
        ],
    )

    history_df = pd.DataFrame(
        history.history
    )

    history_df.to_csv(
        os.path.join(
            output_dir,
            "history.csv",
        ),
        index=False,
    )

    # Restore the model with the highest validation accuracy.
    model.load_weights(
        weight_path
    )

    return (
        model,
        history_df,
    )


# ============================================================
# Confusion matrix
# ============================================================

def make_confusion_matrix(
    y_true,
    y_pred,
    num_classes,
):
    """
    Generate confusion matrix.
    """

    matrix = np.zeros(
        (
            num_classes,
            num_classes,
        ),
        dtype=np.int64,
    )

    for true_index, predicted_index in zip(
        y_true,
        y_pred,
    ):

        matrix[
            true_index,
            predicted_index,
        ] += 1

    return matrix


# ============================================================
# Save confusion matrix
# ============================================================

def save_confusion_matrix(
    matrix,
    trna_names,
    output_dir,
):
    """
    Save raw and normalized confusion matrices
    and corresponding figures.
    """

    raw_df = pd.DataFrame(
        matrix,
        index=trna_names,
        columns=trna_names,
    )

    raw_df.to_csv(
        os.path.join(
            output_dir,
            "confusion_matrix_counts.tsv",
        ),
        sep="\t",
    )

    denominator = matrix.sum(
        axis=1,
        keepdims=True,
    ).astype(
        float
    )

    denominator[
        denominator
        == 0
    ] = 1.0

    normalized = (
        matrix
        / denominator
    )

    normalized_df = pd.DataFrame(
        normalized,
        index=trna_names,
        columns=trna_names,
    )

    normalized_df.to_csv(
        os.path.join(
            output_dir,
            "confusion_matrix_normalized.tsv",
        ),
        sep="\t",
    )

    plt.figure(
        figsize=(
            10,
            9,
        )
    )

    plt.imshow(
        normalized,
        aspect="auto",
        vmin=0.0,
        vmax=1.0,
    )

    plt.colorbar(
        label="Fraction"
    )

    plt.xticks(
        np.arange(
            len(
                trna_names
            )
        ),
        trna_names,
        rotation=90,
    )

    plt.yticks(
        np.arange(
            len(
                trna_names
            )
        ),
        trna_names,
    )

    plt.xlabel(
        "Predicted tRNA"
    )

    plt.ylabel(
        "True tRNA"
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "confusion_matrix.png",
        ),
        dpi=300,
    )

    plt.savefig(
        os.path.join(
            output_dir,
            "confusion_matrix.pdf",
        )
    )

    plt.close()


# ============================================================
# Learning curves
# ============================================================

def plot_learning_curves(
    history_df,
    output_dir,
):
    """
    Save training/validation accuracy and loss curves.
    """

    epochs = (
        np.arange(
            len(
                history_df
            )
        )
        + 1
    )

    # Accuracy curve.
    plt.figure(
        figsize=(
            7,
            5,
        )
    )

    plt.plot(
        epochs,
        history_df[
            "accuracy"
        ],
        label="Training",
        linewidth=1.2,
    )

    plt.plot(
        epochs,
        history_df[
            "val_accuracy"
        ],
        label="Validation",
        linewidth=1.2,
        linestyle="--",
    )

    plt.xlabel(
        "Epoch"
    )

    plt.ylabel(
        "Accuracy"
    )

    plt.xlim(
        1,
        len(
            history_df
        ),
    )

    plt.legend(
        frameon=False
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "learning_curve_accuracy.png",
        ),
        dpi=300,
    )

    plt.savefig(
        os.path.join(
            output_dir,
            "learning_curve_accuracy.pdf",
        )
    )

    plt.close()

    # Loss curve.
    plt.figure(
        figsize=(
            7,
            5,
        )
    )

    plt.plot(
        epochs,
        history_df[
            "loss"
        ],
        label="Training",
        linewidth=1.2,
    )

    plt.plot(
        epochs,
        history_df[
            "val_loss"
        ],
        label="Validation",
        linewidth=1.2,
        linestyle="--",
    )

    plt.xlabel(
        "Epoch"
    )

    plt.ylabel(
        "Loss"
    )

    plt.xlim(
        1,
        len(
            history_df
        ),
    )

    plt.legend(
        frameon=False
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            output_dir,
            "learning_curve_loss.png",
        ),
        dpi=300,
    )

    plt.savefig(
        os.path.join(
            output_dir,
            "learning_curve_loss.pdf",
        )
    )

    plt.close()


# ============================================================
# Evaluate model
# ============================================================

def evaluate_model(
    model_name,
    model,
    history_df,
    x_validation,
    y_validation_index,
    trna_names,
    output_dir,
):
    """
    Evaluate the best checkpoint and generate outputs.
    """

    print()
    print(
        "Evaluating:",
        model_name
    )

    probabilities = model.predict(
        x_validation,
        batch_size=BATCH_SIZE,
        verbose=1,
    )

    y_pred = np.argmax(
        probabilities,
        axis=1,
    )

    matrix = make_confusion_matrix(
        y_validation_index,
        y_pred,
        len(
            trna_names
        ),
    )

    save_confusion_matrix(
        matrix,
        trna_names,
        output_dir,
    )

    plot_learning_curves(
        history_df,
        output_dir,
    )

    validation_accuracy = float(
        np.mean(
            y_pred
            == y_validation_index
        )
    )

    best_epoch = int(
        history_df[
            "val_accuracy"
        ].idxmax()
        + 1
    )

    best_val_accuracy = float(
        history_df[
            "val_accuracy"
        ].max()
    )

    best_val_loss = float(
        history_df.loc[
            best_epoch - 1,
            "val_loss",
        ]
    )

    summary = {
        "model":
            model_name,

        "validation_accuracy_best_weights":
            validation_accuracy,

        "best_epoch":
            best_epoch,

        "best_val_accuracy_history":
            best_val_accuracy,

        "val_loss_at_best_epoch":
            best_val_loss,

        "n_train_per_class":
            N_TRAIN_PER_CLASS,

        "n_validation_per_class":
            N_VALIDATION_PER_CLASS,

        "total_training":
            N_TRAIN_PER_CLASS
            * len(
                trna_names
            ),

        "total_validation":
            N_VALIDATION_PER_CLASS
            * len(
                trna_names
            ),

        "epochs":
            EPOCHS,

        "batch_size":
            BATCH_SIZE,

        "learning_rate":
            LEARNING_RATE,

        "seed":
            SEED,
    }

    pd.DataFrame(
        [
            summary
        ]
    ).to_csv(
        os.path.join(
            output_dir,
            "summary.tsv",
        ),
        sep="\t",
        index=False,
    )

    print()
    print(
        model_name,
        "validation accuracy:",
        validation_accuracy,
    )

    print(
        model_name,
        "best epoch:",
        best_epoch,
    )

    print(
        model_name,
        "best val_accuracy:",
        best_val_accuracy,
    )

    return summary


# ============================================================
# Save model comparison
# ============================================================

def save_comparison(
    summaries,
):
    """
    Save summary table for CNN-only and CNN-WaveNet.
    """

    comparison_df = pd.DataFrame(
        summaries
    )

    output_path = os.path.join(
        OUTPUT_DIR,
        "architecture_comparison.tsv",
    )

    comparison_df.to_csv(
        output_path,
        sep="\t",
        index=False,
    )

    print()
    print(
        "=============================================="
    )

    print(
        "Architecture comparison"
    )

    print(
        "=============================================="
    )

    print()

    print(
        comparison_df[
            [
                "model",
                "validation_accuracy_best_weights",
                "best_epoch",
                "best_val_accuracy_history",
            ]
        ].to_string(
            index=False
        )
    )

    print()
    print(
        "Saved:",
        output_path,
    )


# ============================================================
# Main
# ============================================================

def main():
    """
    Run CNN-only and CNN-WaveNet architecture comparison.
    """

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    print(
        "Visible GPUs:",
        tf.config.list_physical_devices(
            "GPU"
        ),
    )

    set_seed(
        SEED
    )

    (
        x_train,
        y_train,
        x_validation,
        y_validation,
        y_validation_index,
        trna_names,
    ) = load_data(
        INPUT_PARQUET
    )

    save_class_index(
        trna_names
    )

    summaries = []

    for model_name in MODEL_NAMES:

        model_output_dir = os.path.join(
            OUTPUT_DIR,
            model_name,
        )

        os.makedirs(
            model_output_dir,
            exist_ok=True,
        )

        (
            model,
            history_df,
        ) = train_model(
            model_name=model_name,
            x_train=x_train,
            y_train=y_train,
            x_validation=x_validation,
            y_validation=y_validation,
            output_dir=model_output_dir,
        )

        summary = evaluate_model(
            model_name=model_name,
            model=model,
            history_df=history_df,
            x_validation=x_validation,
            y_validation_index=y_validation_index,
            trna_names=trna_names,
            output_dir=model_output_dir,
        )

        summaries.append(
            summary
        )

        # Release memory before training the next model.
        del model

        keras.backend.clear_session()

    save_comparison(
        summaries
    )


if __name__ == "__main__":
    main()