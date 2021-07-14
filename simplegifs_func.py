import imageio

def simplegifs(imglist, savename):
    # let's make a gif!
    giflist = []
    for filename in imglist:
        giflist.append(imageio.imread(filename))
    imageio.mimsave(savename,giflist,duration=1)

"""
imglist = []
for i in range(0,19):
    imglist.append('/Users/rileyannereid/Dropbox/NSET/plots/0720gifs/0521250_pos/' +str(i)+'.png')

savename = '/Users/rileyannereid/Dropbox/NSET/plots/0720gifs/0521250_pos/28kraylocs.gif'
simplegifs(imglist, savename)
"""