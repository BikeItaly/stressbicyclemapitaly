var boundary = [
    [45.151053,8.312531],
    [45.776623,10.062103]
];

const map = L.map('mapid').setView([45.473027,9.173627], 15)
map.setMaxBounds(boundary);
var hash = new L.Hash(map);
const settings = [{ color: '#0099cc', weight: 3, key: 'LTS1', zIndex: 1, title: 'LTS 1 - Adatto ai bambini', url: 'data/level_1.json' },
                { color: '#1C7C54', weight: 3, key: 'LTS2', zIndex: 2, title: 'LTS 2 - Basso Stress', url: 'data/level_2.json' },
                { color: '#F0C808', weight: 3, key: 'LTS3', zIndex: 3, title: 'LTS 3 - Stress Moderato', url: 'data/level_3.json' },
                { color: '#DD5454', weight: 3, key: 'LTS4', zIndex: 4, title: 'LTS 4 - Stress Alto', url: 'data/level_4.json' }]
const homePage = ''
const legendTitle = 'Livello di stress da traffico in bicicletta - Cagliari'
const layers = {}
const tree = rbush.rbush();

addLegend()
addStressLayers()
addIconLayers()


///// Functions ////

function addLegend () {
  const legend = L.control({position: 'topright'})
  legend.onAdd = function (map) {
    const div = L.DomUtil.create('div', 'info legend')
    let legendHtml = '<center><a href="' + homePage + '" target="_blank"><h3>' + legendTitle + '</h3></a></center><table>'
    for (let setting of settings) {
      legendHtml += addLegendLine(setting)
    }
    legendHtml += '</table>'
    div.innerHTML = legendHtml
    div.addEventListener('mouseover', function () {map.doubleClickZoom.disable(); });
    div.addEventListener('mouseout', function () {map.doubleClickZoom.enable(); });
    return div
  }
  legend.addTo(map)
}

function addStressLayers () {
  for (let setting of settings) {
    addStressLayerToMap(setting)
  }
}

function addStressLayerToMap (setting) {
  const xhr = new XMLHttpRequest()
  xhr.open('GET', setting.url)
  xhr.setRequestHeader('Content-Type', 'application/json')
  xhr.onload = function () {
    if (xhr.status === 200) {
      const data = JSON.parse(xhr.responseText)
      const tileIndex = geojsonvt(data, { maxZoom: 18 })
      tree.load(data)

      const canvasTiles = L.tileLayer.canvas()
      canvasTiles.drawTile = function (canvas, tilePoint, zoom) {
        const tile = tileIndex.getTile(zoom, tilePoint.x, tilePoint.y)
        if (!tile) { return }
        drawFeatures(canvas.getContext('2d'), tile.features, setting.color, setting.weight)
      }
      canvasTiles.addTo(map)
      layers[setting.key] = canvasTiles
    } else {
      alert('Request failed.  Returned status of ' + xhr.status)
    }
  }
  xhr.send()
}

function drawFeatures (ctx, features, lineColor, weight) {
  ctx.strokeStyle = lineColor
  ctx.lineWidth = weight

  for (let feature of features) {
    const type = feature.type
    ctx.fillStyle = feature.tags.color ? feature.tags.color : 'rgba(255,0,0,0.05)'
    ctx.beginPath()
    for (let geom of feature.geometry) {
      const pad = 1
      const ratio = .1
      if (type === 1) {
        ctx.arc(geom[0] * ratio + pad, geom[1] * ratio + pad, 2, 0, 2 * Math.PI, false)
        continue
      }
      for (var k = 0; k < geom.length; k++) {
        var p = geom[k]
        var extent = 4096
        var x = p[0] / extent * 256
        var y = p[1] / extent * 256
        if (k) {
          ctx.lineTo(x + pad, y + pad)
        } else {
          ctx.moveTo(x + pad, y + pad)
        }
      }
    }
    if (type === 3 || type === 1) ctx.fill('evenodd')
    ctx.stroke()
  }
}

function toggleLayer (checkbox) {
  if (checkbox.checked) {
    map.addLayer(layers[checkbox.id])
  } else {
    map.removeLayer(layers[checkbox.id])
  }
}

function addLegendLine (setting) {
  return ('<tr><td><input type="checkbox" id="' +
    setting.key +
    '" onclick="toggleLayer(this)" checked /></td>' +
    '<td><hr style="display:inline-block; width: 50px;" color="' +
    setting.color +
    '" size="5" /></td><td>' +
    setting.title +
    '</td></tr>'
  )
}


function addIconLayers(){

  const providers = [];
 /*
  providers.push({
      title: 'JawgMaps',
      icon: 'img/icons-streets.png',
      layer: L.tileLayer('https://tile.jawg.io/{z}/{x}/{y}.png?api-key=community', {
          attribution: "&copy; OpenStreetMap JawgMaps",
          maxZoom: 18
      })
  });
  
  */
  
  providers.push({
      title: 'OpenStreetMap',
      icon: 'img/icons-carto.png',
      layer: L.tileLayer('https://cartodb-basemaps-{s}.global.ssl.fastly.net/dark_all/{z}/{x}/{y}.png', {
          maxZoom: 18,
          attribution: 'Map tiles by [[http://cartodb.com/attributions#basemaps|CartoDB]], under [[https://creativecommons.org/licenses/by/3.0/|CC BY 3.0]]. Data by [[http://www.openstreetmap.org/|OpenStreetMap]], under ODbL.'
      })
  });  
  providers.push({
      title: 'OpenStreetMap',
      icon: 'img/icons-mapnik.png',
      layer: L.tileLayer('http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
          maxZoom: 18,
          attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
      })
  });

  providers.push({
      title: 'OSM b/n',
      icon: 'img/icons-osm-bw.png',
      layer: L.tileLayer('http://{s}.tiles.wmflabs.org/bw-mapnik/{z}/{x}/{y}.png', {
          maxZoom: 18,
          attribution: '&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>'
      })
  });


  providers.push({
      title: 'Fotoaerea Trento',
      icon: 'img/icons-satellite.png',
      layer: L.tileLayer('https://tiles.openaerialmap.org/60770b0fb85cd80007a01414/0/60770b0fb85cd80007a01415/{z}/{x}/{y}', {
          attribution: "&copy; Comune di Trento",
          maxZoom: 18
      })
  });

  L.control.iconLayers(providers).addTo(map);

}


function getFeaturesNearby(point, maxMeters, breakOnFirst)
{
  ret = [];
  const pt = turf.helpers.point(point);
  const nearby = tree.search(pt);
  for(let feature of nearby.features){
    if(breakOnFirst && ret.length){return ret;}
    const line = turf.helpers.lineString(feature.geometry.coordinates);
    if(turf.pointToLineDistance(pt, line, {units: 'meters'})<maxMeters){
      ret.push(feature);
    }
  }

  return ret;
}


function displayOsmElementInfo(element, latlng) {

  const xhr = new XMLHttpRequest()
  xhr.open('GET','https://api.openstreetmap.org/api/0.6/'+element)
  xhr.onload = function () {
    let popup = '<b><a href="https://www.openstreetmap.org/' + element + '" target="_blank">' + element + '</a></b><hr>'
    if (xhr.status === 200) {
      const xmlDOM = new DOMParser().parseFromString(xhr.responseText, 'text/xml');
      const tags = xmlDOM.getElementsByTagName("tag");
      for(let i=0; i<tags.length; i++)
      {
        popup += tags[i].attributes["k"].value+": <b>"+tags[i].attributes["v"].value+'</b><br>';
      }
    } else {
      popup += 'Failed to request details from osm.org';
    }
    map.openPopup(popup, latlng);
  }
  xhr.send()
}


let highlight;
let timer;
map.on('mousemove', function(e) {
  const features = getFeaturesNearby([e.latlng.lng,e.latlng.lat], 5, true)
  clearTimeout(timer);
  if (features.length!=0) {
    document.getElementById('mapid').style.cursor = 'pointer'
  }
  else {
    timer = setTimeout(function()
                {
	                 document.getElementById('mapid').style.cursor = ''
                 }, 100);
  }
})

map.on('click', function(e) {
  if (highlight){
    map.removeLayer(highlight)
  }
  const features = getFeaturesNearby([e.latlng.lng,e.latlng.lat], 5, true);
  if (features.length!=0) {
    displayOsmElementInfo(features[0].id, e.latlng);
    highlight = new L.geoJson(features[0],{style: {color:'#df42f4',  weight: 5}}).addTo(map);
    map.on('popupclose', function() {
     map.removeLayer(highlight)
   });
  }
 });
