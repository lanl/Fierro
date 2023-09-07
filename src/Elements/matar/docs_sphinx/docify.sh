make html 
cp renamer.py _build/html/
cd _build/html
python3 renamer.py
cp -r ./*  ../../../docs/
cp -r ./_static/* ../../../docs/static/
