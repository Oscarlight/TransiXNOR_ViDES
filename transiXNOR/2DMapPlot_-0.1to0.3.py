from skimage.transform import resize
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
matplotlib.use('agg')

# This plot is only valid for Python 3.6. If you use Python 2.7, the color scheme is not correct

model_path             = './D9'

label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
mpl.rcParams['figure.figsize'] = 12, 12

def sqResize(image, outwidth):
    if np.sum(np.abs(image)) > 0:
        log_image = np.log10(np.abs(image))
    else:
        log_image = np.ones_like(image)*(-10)
    abs_image = np.abs(image)
    out_image = resize(log_image, [outwidth, outwidth], mode='constant')
    return out_image

data = np.load(model_path + '/current.npy')
numVds, numVbg, numVtg = data.shape
nRow = 3; nCol = 3; figsize = 8

index = [0, 5, 11, 15, 20, 25, 30, 35, 40]
voltage = [-0.10, -0.05, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
i = 0
for index, vds in zip(index, voltage):
    ax = plt.subplot(nRow, nCol, i+1)
    i += 1
    img = sqResize(data[index, :, :], 41)
    cax = plt.imshow(img, vmin = -2, vmax = 2, origin='lower', interpolation = 'gaussian')
    plt.title(r'$V_{DS}$' + ' = {0:.2f} V'.format(vds), fontsize = label_size)
    plt.ylabel(r'$V_{BG}$ (V)', fontsize = label_size)
    plt.xlabel(r'$V_{TG}$ (V)', fontsize = label_size)
    cbar = plt.colorbar(cax, ticks=[-2.0, -1.0, 0.0, 1.0, 2.0])
    cbar.ax.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$100$'], fontsize = label_size)
    plt.subplots_adjust(wspace = 0.3, hspace = 0.55)
    ax.set_xticklabels([0.0, -0.1, 0.0, 0.1, 0.2, 0.3])
    ax.set_yticklabels([0.0, -0.1, 0.0, 0.1, 0.2, 0.3])

plt.savefig(model_path + '/plots/2Dmapping_-0.1to0.3.pdf')
