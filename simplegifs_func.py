import imageio

# function to turn list of images into a gif

def simplegifs(imglist, savename):
    # let's make a gif!
    giflist = []
    for filename in imglist:
        giflist.append(imageio.imread(filename))
    imageio.mimsave(savename,giflist,duration=1)
