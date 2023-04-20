import sys
import imageio
n = int(sys.argv[1])
dirname = sys.argv[2]
filename = sys.argv[3]
fileout = '{}/{}'.format(dirname, filename)

if n < 366:
    images = []
    for i in range(0, n):
        filein = '{}/{}.png'.format(dirname, i+1)
        images.append(imageio.imread(filein))
    imageio.mimsave(fileout, images)
else:
    with imageio.get_writer(fileout, mode='I') as writer:
        for i in range(0, n):
            filein = '{}/{}.png'.format(dirname, i+1)
            image = imageio.imread(filein)
            writer.append_data(image)
    
