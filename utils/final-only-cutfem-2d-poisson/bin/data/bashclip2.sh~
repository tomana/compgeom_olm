#python data/svgclip.py data/myfile.svg -o data/croppedmyfile.svg
cd "$(dirname "$0")"
cp myfile.svg croppedmyfile.svg

inkscape --verb=EditSelectAll --verb=SelectionUnGroup --verb=EditSelectAll --verb=org.ffaat.filter.fill2stroke.noprefs --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose croppedmyfile.svg
inkscape -f croppedmyfile.svg -A croppedmyfile.pdf
