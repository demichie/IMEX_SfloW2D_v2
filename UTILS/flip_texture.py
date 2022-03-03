import imageio
import numpy as np

im = imageio.imread('texture.png')
print(im.shape)

im_lr = np.fliplr(im)

imageio.imwrite('texture_lr.png', im_lr)

im_ud = np.flipud(im)

imageio.imwrite('texture_ud.png', im_ud)

im_ud_lr = np.flipud(im_lr)

imageio.imwrite('texture_ud_lr.png', im_ud_lr)

print(im.shape)

im_T = np.transpose(im, (1, 0, 2))

imageio.imwrite('texture_T.png', im_T)

im_lr = np.fliplr(im_T)

imageio.imwrite('textureT_lr.png', im_lr)

im_ud = np.flipud(im_T)

imageio.imwrite('textureT_ud.png', im_ud)

im_ud_lr = np.flipud(im_lr)

imageio.imwrite('textureT_ud_lr.png', im_ud_lr)


