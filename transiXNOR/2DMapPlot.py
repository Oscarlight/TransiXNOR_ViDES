import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from skimage.transform import resize

label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 
mpl.rcParams['figure.figsize'] = 12, 12

model_path = './D6'

def sqResize(image, outwidth):
    try:
        log_image = np.log10(np.abs(image))
    except:
        log_image = np.zeros_like(image)
    abs_image = np.abs(image)
    out_image = resize(log_image, [outwidth, outwidth])
    return out_image

data = np.load(model_path + '/current.npy')
numVds, numVbg, numVtg = data.shape
nRow = 3; nCol = 3; figsize = 8
fontsize = label_size

index = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 19, 20]
voltage = [0.00, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.19, 0.20]
i = 0
for index, vds in zip(index, voltage):
    ax = plt.subplot(4, 3, i+1)
    i += 1
    img = sqResize(data[index, :, :], 41)
    cax = plt.imshow(img, vmin = -3, vmax = 2, origin='lower', interpolation = 'gaussian')
    plt.title(r'$V_{DS}$' + ' = {0:.2f} V'.format(vds), fontsize = label_size)
    plt.ylabel(r'$V_{BG}$ (V)', fontsize = label_size)
    plt.xlabel(r'$V_{TG}$ (V)', fontsize = label_size)
    cbar = plt.colorbar(cax, ticks=[-2.0, -1.0, 0.0, 1.0, 2.0])
    cbar.ax.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10$', r'$10^{2}$'], fontsize = label_size)
    plt.subplots_adjust(wspace = 0.3, hspace = 0.55)
    ax.set_xticklabels([0.0, 0.0, 0.1, 0.2])
    ax.set_yticklabels([0.0, 0.0, 0.1, 0.2])

plt.savefig(model_path+'/plots/2DMapping.pdf')
