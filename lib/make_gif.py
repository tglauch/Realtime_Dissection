import imageio
import os
images = []
basepath = '/scratch9/tglauch/realtime_service/output/IC190504/4FGL_J0420.3-3745/lightcurve/'
filenames = sorted([os.path.join(basepath, i, 'sed.pdf') for i in os.listdir(basepath) if
                    os.path.isdir(os.path.join(basepath, i))])
print filenames
for filename in filenames:
    os.system('convert -density 150 {} -quality 90 {}'.format(filename, filename.replace('.pdf', '.png')))
    images.append(imageio.imread(filename.replace('.pdf', '.png')))
imageio.mimsave('./movie.gif', images)
