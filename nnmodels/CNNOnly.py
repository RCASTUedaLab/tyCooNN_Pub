import tensorflow as tf

from tensorflow.keras.layers import (
    Input,
    Conv1D,
    MaxPooling1D,
    AveragePooling1D,
    BatchNormalization,
    Dropout,
    Activation,
    Concatenate,
    Add,
    GaussianNoise,
    GlobalAveragePooling1D,
    Softmax,
)
from tensorflow.keras.models import Model
from tensorflow.keras.regularizers import l2

from funcy import identity, juxt, rcompose


# ==========================================================
# Mish activation
# Compatible with older TensorFlow / Keras
# ==========================================================

def mish(x):
    return x * tf.math.tanh(
        tf.math.softplus(x)
    )


# ==========================================================
# Basic convolution layers
# ==========================================================

def conv1D(filters, kernel_size):
    return Conv1D(
        filters,
        kernel_size,
        padding="same",
        kernel_initializer="he_normal",
        kernel_regularizer=l2(0.0001),
    )


def conv1D_halve(filters, kernel_size):
    return Conv1D(
        filters,
        kernel_size,
        padding="same",
        strides=2,
        kernel_initializer="he_normal",
        kernel_regularizer=l2(0.0001),
    )


def convBlock(
    f1,
    k1,
    f2,
    k2,
    f3,
    k3,
    do_r,
):
    return rcompose(
        conv1D(f1, k1),
        conv1D(f2, k2),
        conv1D(f3, k3),
        MaxPooling1D(pool_size=2),
        BatchNormalization(),
        Dropout(do_r),
    )


# ==========================================================
# CNN feature extractor
# This corresponds to the CNN part of the current model
# before the WaveNet layers.
# ==========================================================

def build_cnn_feature_extractor(
    inputs,
    do_r=0.2,
):

    def ljuxt(*fs):
        return rcompose(
            juxt(*fs),
            list,
        )

    def fire_module(
        filters_squeeze,
        filters_expand,
    ):
        return rcompose(
            BatchNormalization(),
            Activation(mish),
            conv1D(
                filters_squeeze,
                1,
            ),
            BatchNormalization(),
            Activation(mish),
            ljuxt(
                conv1D(
                    filters_expand // 2,
                    1,
                ),
                conv1D(
                    filters_expand // 2,
                    3,
                ),
            ),
            Concatenate(),
        )

    def fire_module_with_shortcut(
        filters_squeeze,
        filters_expand,
    ):
        return rcompose(
            ljuxt(
                fire_module(
                    filters_squeeze,
                    filters_expand,
                ),
                identity,
            ),
            Add(),
        )

    def inception():

        u1 = rcompose(
            AveragePooling1D(
                pool_size=3,
                strides=1,
                padding="same",
            ),
            conv1D(
                48,
                1,
            ),
        )

        u2 = conv1D(
            48,
            1,
        )

        u3 = rcompose(
            conv1D(
                16,
                1,
            ),
            conv1D(
                48,
                3,
            ),
        )

        u4 = rcompose(
            conv1D(
                16,
                1,
            ),
            conv1D(
                48,
                3,
            ),
            conv1D(
                48,
                3,
            ),
        )

        return rcompose(
            ljuxt(
                u1,
                u2,
                u3,
                u4,
            ),
            Concatenate(axis=2),
        )

    cnn_block = rcompose(
        GaussianNoise(
            stddev=0.02
        ),

        conv1D_halve(
            48,
            3,
        ),

        BatchNormalization(),

        Dropout(
            do_r
        ),

        convBlock(
            48, 3,
            48, 3,
            48, 3,
            do_r,
        ),

        convBlock(
            16, 1,
            48, 3,
            48, 3,
            do_r,
        ),

        convBlock(
            16, 1,
            48, 3,
            48, 3,
            do_r,
        ),

        fire_module(
            8,
            64,
        ),

        fire_module_with_shortcut(
            8,
            64,
        ),

        fire_module(
            16,
            128,
        ),

        MaxPooling1D(),

        convBlock(
            48, 3,
            48, 3,
            48, 3,
            do_r,
        ),

        fire_module(
            8,
            64,
        ),

        fire_module_with_shortcut(
            8,
            64,
        ),

        fire_module(
            16,
            128,
        ),

        MaxPooling1D(),

        fire_module_with_shortcut(
            16,
            128,
        ),

        BatchNormalization(),

        Dropout(
            do_r
        ),

        Activation(
            mish
        ),

        inception(),
    )

    return cnn_block(
        inputs
    )


# ==========================================================
# CNN-only classifier
# ==========================================================

def build_network(
    shape,
    num_classes,
    do_r=0.2,
):

    inputs = Input(
        batch_shape=shape
    )

    x = build_cnn_feature_extractor(
        inputs,
        do_r=do_r,
    )

    x = Conv1D(
        num_classes,
        1,
        padding="same",
    )(
        x
    )

    x = GlobalAveragePooling1D()(
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