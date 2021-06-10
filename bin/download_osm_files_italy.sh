url="https://download.geofabrik.de/europe/italy/"
nord_est="nord-est-latest.osm.pbf"
nord_ovest="nord-ovest-latest.osm.pbf"
sud="sud-latest.osm.pbf"
isole="isole-latest.osm.pbf"
centro="centro-latest.osm.pbf"
wget $url$nord_est
sourcefile=$nord_est
osmium tags-filter --overwrite "${sourcefile}" w/barrier -o barriers.pbf
osmium tags-filter --overwrite "${sourcefile}" w/highway -o highways.pbf
osmium tags-filter --overwrite -i highways.pbf 'w/footway=sidewalk' -o highwaysnotsidewalk.pbf
osmium tags-filter --overwrite "${sourcefile}" w/footway=sidewalk -o footways.pbf
osmium tags-filter --overwrite footways.pbf w/bicycle -o footways_bicycle.pbf
osmium tags-filter --overwrite -i footways_bicycle.pbf w/bicycle=no -o footways_bicycle_yes.pbf
osmium merge --overwrite barriers.pbf highwaysnotsidewalk.pbf footways_bicycle_yes.pbf -o bike_nord_est.osm
wget $url$nord_ovest
sourcefile=$nord_ovest
osmium tags-filter --overwrite "${sourcefile}" w/barrier -o barriers.pbf
osmium tags-filter --overwrite "${sourcefile}" w/highway -o highways.pbf
osmium tags-filter --overwrite -i highways.pbf 'w/footway=sidewalk' -o highwaysnotsidewalk.pbf
osmium tags-filter --overwrite "${sourcefile}" w/footway=sidewalk -o footways.pbf
osmium tags-filter --overwrite footways.pbf w/bicycle -o footways_bicycle.pbf
osmium tags-filter --overwrite -i footways_bicycle.pbf w/bicycle=no -o footways_bicycle_yes.pbf
osmium merge --overwrite barriers.pbf highwaysnotsidewalk.pbf footways_bicycle_yes.pbf -o bike_nord_ovest.osmwget $url$nord_ovest
wget $url$sud
sourcefile=$sud
osmium tags-filter --overwrite "${sourcefile}" w/barrier -o barriers.pbf
osmium tags-filter --overwrite "${sourcefile}" w/highway -o highways.pbf
osmium tags-filter --overwrite -i highways.pbf 'w/footway=sidewalk' -o highwaysnotsidewalk.pbf
osmium tags-filter --overwrite "${sourcefile}" w/footway=sidewalk -o footways.pbf
osmium tags-filter --overwrite footways.pbf w/bicycle -o footways_bicycle.pbf
osmium tags-filter --overwrite -i footways_bicycle.pbf w/bicycle=no -o footways_bicycle_yes.pbf
osmium merge --overwrite barriers.pbf highwaysnotsidewalk.pbf footways_bicycle_yes.pbf -o bike_sud.osm
wget $url$isole
sourcefile=$isole
osmium tags-filter --overwrite "${sourcefile}" w/barrier -o barriers.pbf
osmium tags-filter --overwrite "${sourcefile}" w/highway -o highways.pbf
osmium tags-filter --overwrite -i highways.pbf 'w/footway=sidewalk' -o highwaysnotsidewalk.pbf
osmium tags-filter --overwrite "${sourcefile}" w/footway=sidewalk -o footways.pbf
osmium tags-filter --overwrite footways.pbf w/bicycle -o footways_bicycle.pbf
osmium tags-filter --overwrite -i footways_bicycle.pbf w/bicycle=no -o footways_bicycle_yes.pbf
osmium merge --overwrite barriers.pbf highwaysnotsidewalk.pbf footways_bicycle_yes.pbf -o bike_isole.osm
wget $url$centro
sourcefile=$centro
osmium tags-filter --overwrite "${sourcefile}" w/barrier -o barriers.pbf
osmium tags-filter --overwrite "${sourcefile}" w/highway -o highways.pbf
osmium tags-filter --overwrite -i highways.pbf 'w/footway=sidewalk' -o highwaysnotsidewalk.pbf
osmium tags-filter --overwrite "${sourcefile}" w/footway=sidewalk -o footways.pbf
osmium tags-filter --overwrite footways.pbf w/bicycle -o footways_bicycle.pbf
osmium tags-filter --overwrite -i footways_bicycle.pbf w/bicycle=no -o footways_bicycle_yes.pbf
osmium merge --overwrite barriers.pbf highwaysnotsidewalk.pbf footways_bicycle_yes.pbf -o bike_centro.osm