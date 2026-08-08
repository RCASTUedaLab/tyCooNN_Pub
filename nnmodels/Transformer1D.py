import tensorflow as tf

from tensorflow.keras.layers import (
    Input,
    Conv1D,
    Dense,
    Dropout,
    LayerNormalization,
    MultiHeadAttention,
    GlobalAveragePooling1D,
    Softmax,
)
from tensorflow.keras.models import Model
from tensorflow.keras.regularizers import l2

from nnmodels.CNNOnly import (
    build_cnn_feature_extractor,
    mish,
)


# ==========================================================
# Positional embedding
# ==========================================================

class TrainablePositionEmbedding(
    tf.keras.layers.Layer
):

    def build(
        self,
        input_shape,
    ):

        sequence_length = int(
            input_shape[1]
        )

        embedding_dim = int(
            input_shape[2]
        )

        self.position_embedding = (
            self.add_weight(
                name="position_embedding",
                shape=(
                    sequence_length,
                    embedding_dim,
                ),
                initializer="random_normal",
                trainable=True,
            )
        )

        super(
            TrainablePositionEmbedding,
            self
        ).build(
            input_shape
        )

    def call(
        self,
        inputs,
    ):

        return (
            inputs
            + self.position_embedding
        )


# ==========================================================
# Transformer encoder
# ==========================================================

def transformer_block(
    x,
    embed_dim,
    num_heads,
    ff_dim,
    dropout_rate,
):

    # ------------------------------------------------------
    # Multi-head self-attention
    # ------------------------------------------------------

    attention_output = MultiHeadAttention(
        num_heads=num_heads,
        key_dim=embed_dim // num_heads,
        dropout=dropout_rate,
    )(
        x,
        x,
    )

    attention_output = Dropout(
        dropout_rate
    )(
        attention_output
    )

    x = LayerNormalization(
        epsilon=1e-6
    )(
        x
        + attention_output
    )

    # ------------------------------------------------------
    # Feed-forward network
    # ------------------------------------------------------

    ff_output = Dense(
        ff_dim,
        activation=mish,
        kernel_regularizer=l2(
            0.0001
        ),
    )(
        x
    )

    ff_output = Dropout(
        dropout_rate
    )(
        ff_output
    )

    ff_output = Dense(
        embed_dim,
        kernel_regularizer=l2(
            0.0001
        ),
    )(
        ff_output
    )

    ff_output = Dropout(
        dropout_rate
    )(
        ff_output
    )

    x = LayerNormalization(
        epsilon=1e-6
    )(
        x
        + ff_output
    )

    return x


# ==========================================================
# Transformer classifier
# ==========================================================

def build_network(
    shape,
    num_classes,
    do_r=0.2,
    embed_dim=128,
    num_heads=4,
    ff_dim=256,
    num_transformer_blocks=2,
):

    inputs = Input(
        batch_shape=shape
    )

    # ------------------------------------------------------
    # Use exactly the same CNN feature extractor
    # as the CNN-only and CNN-WaveNet architectures.
    # ------------------------------------------------------

    x = build_cnn_feature_extractor(
        inputs,
        do_r=do_r,
    )

    # ------------------------------------------------------
    # Project CNN features to Transformer embedding size.
    # ------------------------------------------------------

    x = Conv1D(
        embed_dim,
        1,
        padding="same",
        kernel_regularizer=l2(
            0.0001
        ),
    )(
        x
    )

    # ------------------------------------------------------
    # Positional information
    # ------------------------------------------------------

    x = TrainablePositionEmbedding()(
        x
    )

    # ------------------------------------------------------
    # Transformer encoders
    # ------------------------------------------------------

    for _ in range(
        num_transformer_blocks
    ):

        x = transformer_block(
            x=x,
            embed_dim=embed_dim,
            num_heads=num_heads,
            ff_dim=ff_dim,
            dropout_rate=do_r,
        )

    # ------------------------------------------------------
    # Classification
    # ------------------------------------------------------

    x = GlobalAveragePooling1D()(
        x
    )

    x = Dense(
        128,
        activation=mish,
        kernel_regularizer=l2(
            0.0001
        ),
    )(
        x
    )

    x = Dropout(
        do_r
    )(
        x
    )

    x = Dense(
        num_classes
    )(
        x
    )

    outputs = Softmax()(
        x
    )

    model = Model(
        inputs=inputs,
        outputs=outputs,
    )

    return model