#!/bin/bash
#echo "level_1"
#cp centro/level_1.json geojson/output.geojson
cd isole
#ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_1.json -nln output
cd ../nord-est
#ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_1.json -nln output
cd ../nord-ovest
#ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_1.json -nln output
cd ../sud
#ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_1.json -nln output
cd ../geojson
#tippecanoe -e level_1 -S10 -al -an -ad -aC -g2 -aG -pn -pC -pS -z16 -Z6 output.geojson
#mv output.geojson level_1_italy.geojson
#zip level_1_italy.zip level_1_italy.geojson
#rm level_1_italy.geojson
echo "level_2"
cd ..
cp centro/level_2.json geojson/output.geojson
cd isole
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_2.json -nln output
cd ../nord-est
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_2.json -nln output
cd ../nord-ovest
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_2.json -nln output
cd ../sud
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_2.json -nln output
cd ../geojson
#tippecanoe -e level_2 -S10 -al -an -ad -aC -g2 -aG -pn -pC -pS -z16 -Z6 output.geojson
mv output.geojson level_2_italy.geojson
#zip level_2_italy.zip level_2_italy.geojson
#rm level_2_italy.geojson
echo "level_3"
cd ..
cp centro/level_3.json geojson/output.geojson
cd isole
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_3.json -nln output
cd ../nord-est
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_3.json -nln output
cd ../nord-ovest
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_3.json -nln output
cd ../sud
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_3.json -nln output
cd ../geojson
#tippecanoe -e level_3 -S10 -al -an -ad -aC -g2 -aG -pn -pC -pS -z16 -Z6 output.geojson
mv output.geojson level_3_italy.geojson
#zip level_3_italy.zip level_3_italy.geojson
#rm level_3_italy.geojson
echo "level_4"
cd ..
cp centro/level_4.json geojson/output.geojson
cd isole
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_4.json -nln output
cd ../nord-est
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_4.json -nln output
cd ../nord-ovest
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_4.json -nln output
cd ../sud
ogr2ogr -f 'GeoJSON' -update -append ../geojson/output.geojson level_4.json -nln output
cd ../geojson
#tippecanoe -e level_4 -S10 -al -an -ad -aC -g2 -aG -pn -pC -pS -z16 -Z6 output.geojson
mv output.geojson level_4_italy.geojson
#zip level_4_italy.zip level_4_italy.geojson
#rm level_4_italy.geojson
