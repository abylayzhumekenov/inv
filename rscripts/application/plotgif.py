import imageio
images = []
for i in range(0, 365):
    filename = 'img/{}.png'.format(i)
    images.append(imageio.imread(filename))
imageio.mimsave('field.gif', images)
