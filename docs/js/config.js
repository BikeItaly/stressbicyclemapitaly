var MAPBOX_TOKEN = 'pk.eyJ1IjoibmFwbyIsImEiOiJ5dG5HeDFZIn0.v4diQq5k3tVzeyRdjYks_A';
var MAPBOX_STYLE = 'mapbox://styles/mapbox/light-v10';
var STRESSCYCLELEVES_TILESET = 'mapbox://napo.stress_cycle_levels'
var STRESSCYCLELEVES_SOURCE = 'stress_cycle_levels'
var STRESSCYCLELEVES_LAYER = 'stresslevels'
var STRESSCYCLELEVES_SOURCEmin = 'stresslevels_layer'
var STRESSCYCLELEVES_TILESETmin = 'https://bicistressatedaltraffico.it/stresslevels_layer/{z}/{x}/{y}.pbf'
var STRESSCYCLELEVES_LAYERmin ='stress_traffic_roads_for_bicycles'
LEVELSANDCOLORS = ['level_1', '#0099cc', 'level_2', '#1C7C54','level_3','#F0C808','level_4','#DD5454'];


var levels={};
levels["level_1"]={};
levels["level_1"].visible=1;
levels["level_1"].color='#0099cc';
levels["level_1"].label="LTS 1 - Adatto ai bambini"
levels["level_2"]={};
levels["level_2"].visible=1;
levels["level_2"].color='#1C7C54';
levels["level_2"].label="LTS 2 - Basso Stress"
levels["level_3"]={};
levels["level_3"].visible=1;
levels["level_3"].color='#F0C808';
levels["level_3"].label="LTS 3 - Stress Moderato"
levels["level_4"]={};
levels["level_4"].visible=1;
levels["level_4"].color='#DD5454';
levels["level_4"].label="LTS 4 - Stress Alto"
