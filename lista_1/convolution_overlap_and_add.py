import numpy as np
from scipy.linalg import toeplitz


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def conv_toeplitz(signal_1, length_signal_2):
    """Turns signal_1 into a Toeplitz matrix X to so you can compute the
    linear convolution signal * signal_2 as np.dot(X, signal_2).

    Args:
        signal_1 (array): Vector that is turned into a Toeplitz matrix used for
        convolution.
        signal_2 (array): The other member of the convolution operation.

    Returns:
        matrix: signal_1 in Toeplitz matrix form.
    """

    M = len(signal_1)
    num_out = length_signal_2 + M - 1

    rowToe = np.append(signal_1[0], np.zeros((1, num_out - M)))
    colToe = np.append(signal_1, np.zeros((num_out - M, 1)))
    return toeplitz(colToe, rowToe)


def overlap_and_add(signal_1, signal_2):

    # x is always the bigger vector
    if len(signal_1) > len(signal_2):
        x = signal_1
        h = signal_2
    else:
        x = signal_2
        h = signal_1

    tam_x = len(x)
    tam_bloco = len(h)
    tam_result = tam_x + tam_bloco - 1

    qtd_bloco = tam_x / tam_bloco

    lista_blocos = list(chunks(x, tam_bloco))

    conv_blocos = []
    for i in range(int(qtd_bloco)):
        conv_matrix = conv_toeplitz(lista_blocos[i], len(h))
        conv_result = np.dot(conv_matrix, h)
        conv_result_with_zeros = np.append(conv_result, np.zeros(tam_x))
        conv_blocos += [conv_result_with_zeros]

    all_conv_blocks = np.hstack(conv_blocos)
    lista_blocks_all = list(chunks(all_conv_blocks, tam_result))
    lista_blocks_only_usable = lista_blocks_all[:(int(np.floor(len(all_conv_blocks) / tam_result)))]
    res = np.sum(lista_blocks_only_usable, axis=0)

    return res


signal_2 = [1, 2, 3, 4, 5, 7]
signal_1 = [1, 2]

resp = overlap_and_add(signal_1, signal_2)
print(resp)
resp_numpy = np.convolve(signal_1, signal_2)
print(resp_numpy)
print(np.all(resp == resp_numpy))
