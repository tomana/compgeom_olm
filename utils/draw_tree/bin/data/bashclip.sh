#python data/svgclip.py data/myfile.svg -o data/croppedmyfile.svg
cd "$(dirname "$0")"
#cp myfile.svg croppedmyfile.svg

python svgclip.py myfile.svg -o croppedmyfile.svg
#rsvg-convert -f pdf -o croppedmyfile.pdf croppedmyfile.svg
inkscape -f croppedmyfile.svg -A croppedmyfile.pdf
