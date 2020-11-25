import tensorflow as tf
import numpy as np

weight_init = tf.contrib.layers.variance_scaling_initializer()
weight_regularizer = tf.contrib.layers.l2_regularizer(0.00001)

def get_model(img_size, CH0, num_classes):
    X    = tf.placeholder(tf.float32, [None, img_size, img_size, 3], name='X')
    Y    = tf.placeholder(tf.float32, [None, num_classes], name='Y')
    LOSS = tf.placeholder(tf.float32, (), name='LOSS')
    KEEP_PROB = tf.placeholder_with_default(1., (), 'KEEP_PROB')
    IS_TRAINING = tf.placeholder(tf.bool, name='IS_TRAINING')
    LR = tf.placeholder(tf.float32, name='learning_rate')

    ch = CH0
    x = X

    L11 = layer(x, CH0*1, is_training=IS_TRAINING, scope=str(1) + '_' + str(1))
    L12 = layer(L11, CH0*1, is_training=IS_TRAINING, scope=str(1) + '_' + str(2))
    L2  = layer(L12, CH0*2, stride=2, is_training=IS_TRAINING, scope=str(2))

    L31 = layer(L2, CH0*2, is_training=IS_TRAINING, scope=str(3) + '_' + str(1))
    L32 = layer(L31, CH0*2, is_training=IS_TRAINING, scope=str(3) + '_' + str(2))
    L4  = layer(L32, CH0*4, stride=2, is_training=IS_TRAINING, scope=str(4))

    L5  = layer(L4, CH0*4, is_training=IS_TRAINING, scope='5')
    L6  = layer(L5, CH0*4, kernel=1, is_training=IS_TRAINING, scope='6')

    x = tf.reduce_mean(L6, axis=[1, 2], keepdims=True)

    x_shape = x.shape.as_list()


    x_rs = tf.reshape(x, [-1, x_shape[-1]])
            
    logits = tf.layers.dense(x_rs, num_classes)
    pred = tf.nn.softmax(logits)

    Loss = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=logits, labels=Y))
    Reg_loss = tf.losses.get_regularization_loss()

    opt = tf.train.MomentumOptimizer(LR, momentum=0.9)
    train_opt = opt.minimize(Loss + Reg_loss)
    
    layers = [L11, L12, L2, L31, L32, L4, L5, L6, x_rs]
    
    return layers, pred, X, Y, IS_TRAINING, LR, train_opt, Loss

def wreg(w):
    return tf.reduce_sum(tf.square(w))


def to_one_hot(array, N):
    array = np.array(array).astype(int)
    labels = np.zeros((array.shape[0], N), dtype=np.float32)
    labels[np.arange(array.shape[0]), array] = 1.
    return labels


def batch_norm(x, is_training=True, scope='batch_norm'):
    return tf.contrib.layers.batch_norm(x,
                                        decay=0.9, epsilon=1e-05,
                                        center=True, scale=True, updates_collections=None,
                                        is_training=is_training, scope=scope)


def conv(x, channels, kernel, stride, padding='SAME', use_bias=True, scope='conv_0'):
    with tf.variable_scope(scope):
        x = tf.layers.conv2d(inputs=x, filters=channels,
                             kernel_size=kernel, kernel_initializer=weight_init,
                             kernel_regularizer=weight_regularizer,
                             strides=stride, use_bias=use_bias, padding=padding)
        return x

    
def layer(x, channels, kernel=3, stride=1, is_training=True, scope='layer'):
    with tf.variable_scope(scope) :

        x = conv(x, channels, kernel=kernel, stride=stride, scope='conv')
        x = batch_norm(x, is_training, scope='batch_norm')
        x = tf.nn.relu(x)
        return x
        
def random_crop(batch, crop_shape=[32,32], padding=4):
    oshape = np.shape(batch[0])

    if padding:
        oshape = (oshape[0] + 2 * padding, oshape[1] + 2 * padding)
    new_batch = []
    npad = ((padding, padding), (padding, padding), (0, 0))
    for i in range(len(batch)):
        new_batch.append(batch[i])
        if padding:
            new_batch[i] = np.lib.pad(batch[i], pad_width=npad,
                                      mode='constant', constant_values=0)
        nh = np.random.randint(0, oshape[0] - crop_shape[0])
        nw = np.random.randint(0, oshape[1] - crop_shape[1])
        new_batch[i] = new_batch[i][nh:nh + crop_shape[0],
                       nw:nw + crop_shape[1]]
    return np.array(new_batch)
