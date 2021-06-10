var MAPBOX_TOKEN = 'pk.eyJ1IjoibmFwbyIsImEiOiJ5dG5HeDFZIn0.v4diQq5k3tVzeyRdjYks_A';
var MAPBOX_STYLE = 'mapbox://styles/mapbox/streets-v11';
var STRESSCYCLELEVES_TILESET = 'mapbox://napo.stress_cycle_levels'
var STRESSCYCLELEVES_SOURCE = 'stress_cycle_levels'
var STRESSCYCLELEVES_LAYER = 'stresslevels'
LEVELSANDCOLORS = ['level_1', '#0099cc', 'level_2', '#1C7C54','level_3','#F0C808','level_4','#DD5454'];

var levels={};
levels["level_1"]={};
levels["level_1"].visible=1;
levels["level_1"].color='#0099cc';
levels["level_2"]={};
levels["level_2"].visible=1;
levels["level_2"].color='#1C7C54';
levels["level_3"]={};
levels["level_3"].visible=1;
levels["level_3"].color='#F0C808';
levels["level_4"]={};
levels["level_4"].visible=1;
levels["level_4"].color='#DD5454';
console.log(levels);
