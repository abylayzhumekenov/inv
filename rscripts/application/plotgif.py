import sys
import imageio
n = int(sys.argv[1])

if n < 365:
    images = []
    for i in range(0, n):
        filename = 'img/{}.png'.format(i+1)
        images.append(imageio.imread(filename))
    imageio.mimsave('field.gif', images)
else:
    with imageio.get_writer('field.gif', mode='I') as writer:
        for i in range(0, n):
            filename = 'img/{}.png'.format(i+1)
            image = imageio.imread(filename)
            writer.append_data(image)
    
