(function(f){if(typeof exports==="object"&&typeof module!=="undefined"){module.exports=f()}else if(typeof define==="function"&&define.amd){define([],f)}else{var g;if(typeof window!=="undefined"){g=window}else if(typeof global!=="undefined"){g=global}else if(typeof self!=="undefined"){g=self}else{g=this}g.geojsonvt = f()}})(function(){var define,module,exports;return (function e(t,n,r){function s(o,u){if(!n[o]){if(!t[o]){var a=typeof require=="function"&&require;if(!u&&a)return a(o,!0);if(i)return i(o,!0);var f=new Error("Cannot find module '"+o+"'");throw f.code="MODULE_NOT_FOUND",f}var l=n[o]={exports:{}};t[o][0].call(l.exports,function(e){var n=t[o][1][e];return s(n?n:e)},l,l.exports,e,t,n,r)}return n[o].exports}var i=typeof require=="function"&&require;for(var o=0;o<r.length;o++)s(r[o]);return s})({1:[function(require,module,exports){
'use strict';

module.exports = clip;

/* clip features between two axis-parallel lines:
 *     |        |
 *  ___|___     |     /
 * /   |   \____|____/
 *     |        |
 */

function clip(features, scale, k1, k2, axis, intersect, minAll, maxAll) {

    k1 /= scale;
    k2 /= scale;

    if (minAll >= k1 && maxAll <= k2) return features; // trivial accept
    else if (minAll > k2 || maxAll < k1) return null; // trivial reject

    var clipped = [];

    for (var i = 0; i < features.length; i++) {

        var feature = features[i],
            geometry = feature.geometry,
            type = feature.type,
            min, max;

        min = feature.min[axis];
        max = feature.max[axis];

        if (min >= k1 && max <= k2) { // trivial accept
            clipped.push(feature);
            continue;
        } else if (min > k2 || max < k1) continue; // trivial reject

        var slices = type === 1 ?
                clipPoints(geometry, k1, k2, axis) :
                clipGeometry(geometry, k1, k2, axis, intersect, type === 3);

        if (slices.length) {
            // if a feature got clipped, it will likely get clipped on the next zoom level as well,
            // so there's no need to recalculate bboxes
            clipped.push({
                geometry: slices,
                type: type,
                tags: features[i].tags || null,
                min: feature.min,
                max: feature.max
            });
        }
    }

    return clipped.length ? clipped : null;
}

function clipPoints(geometry, k1, k2, axis) {
    var slice = [];

    for (var i = 0; i < geometry.length; i++) {
        var a = geometry[i],
            ak = a[axis];

        if (ak >= k1 && ak <= k2) slice.push(a);
    }
    return slice;
}

function clipGeometry(geometry, k1, k2, axis, intersect, closed) {

    var slices = [];

    for (var i = 0; i < geometry.length; i++) {

        var ak = 0,
            bk = 0,
            b = null,
            points = geometry[i],
            area = points.area,
            dist = points.dist,
            len = points.length,
            a, j, last;

        var slice = [];

        for (j = 0; j < len - 1; j++) {
            a = b || points[j];
            b = points[j + 1];
            ak = bk || a[axis];
            bk = b[axis];

            if (ak < k1) {

                if ((bk > k2)) { // ---|-----|-->
                    slice.push(intersect(a, b, k1), intersect(a, b, k2));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk >= k1) slice.push(intersect(a, b, k1)); // ---|-->  |

            } else if (ak > k2) {

                if ((bk < k1)) { // <--|-----|---
                    slice.push(intersect(a, b, k2), intersect(a, b, k1));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk <= k2) slice.push(intersect(a, b, k2)); // |  <--|---

            } else {

                slice.push(a);

                if (bk < k1) { // <--|---  |
                    slice.push(intersect(a, b, k1));
                    if (!closed) slice = newSlice(slices, slice, area, dist);

                } else if (bk > k2) { // |  ---|-->
                    slice.push(intersect(a, b, k2));
                    if (!closed) slice = newSlice(slices, slice, area, dist);
                }
                // | --> |
            }
        }

        // add the last point
        a = points[len - 1];
        ak = a[axis];
        if (ak >= k1 && ak <= k2) slice.push(a);

        // close the polygon if its endpoints are not the same after clipping

        last = slice[slice.length - 1];
        if (closed && last && (slice[0][0] !== last[0] || slice[0][1] !== last[1])) slice.push(slice[0]);

        // add the final slice
        newSlice(slices, slice, area, dist);
    }

    return slices;
}

function newSlice(slices, slice, area, dist) {
    if (slice.length) {
        // we don't recalculate the area/length of the unclipped geometry because the case where it goes
        // below the visibility threshold as a result of clipping is rare, so we avoid doing unnecessary work
        slice.area = area;
        slice.dist = dist;

        slices.push(slice);
    }
    return [];
}

},{}],2:[function(require,module,exports){
'use strict';

module.exports = convert;

var simplify = require('./simplify');

// converts GeoJSON feature into an intermediate projected JSON vector format with simplification data

function convert(data, tolerance) {
    var features = [];

    if (data.type === 'FeatureCollection') {
        for (var i = 0; i < data.features.length; i++) {
            convertFeature(features, data.features[i], tolerance);
        }
    } else if (data.type === 'Feature') {
        convertFeature(features, data, tolerance);

    } else {
        // single geometry or a geometry collection
        convertFeature(features, {geometry: data}, tolerance);
    }
    return features;
}

function convertFeature(features, feature, tolerance) {
    if (feature.geometry === null) {
        // ignore features with null geometry
        return;
    }

    var geom = feature.geometry,
        type = geom.type,
        coords = geom.coordinates,
        tags = feature.properties,
        i, j, rings;

    if (type === 'Point') {
        features.push(create(tags, 1, [projectPoint(coords)]));

    } else if (type === 'MultiPoint') {
        features.push(create(tags, 1, project(coords)));

    } else if (type === 'LineString') {
        features.push(create(tags, 2, [project(coords, tolerance)]));

    } else if (type === 'MultiLineString' || type === 'Polygon') {
        rings = [];
        for (i = 0; i < coords.length; i++) {
            rings.push(project(coords[i], tolerance));
        }
        features.push(create(tags, type === 'Polygon' ? 3 : 2, rings));

    } else if (type === 'MultiPolygon') {
        rings = [];
        for (i = 0; i < coords.length; i++) {
            for (j = 0; j < coords[i].length; j++) {
                rings.push(project(coords[i][j], tolerance));
            }
        }
        features.push(create(tags, 3, rings));

    } else if (type === 'GeometryCollection') {
        for (i = 0; i < geom.geometries.length; i++) {
            convertFeature(features, {
                geometry: geom.geometries[i],
                properties: tags
            }, tolerance);
        }

    } else {
        throw new Error('Input data is not a valid GeoJSON object.');
    }
}

function create(tags, type, geometry) {
    var feature = {
        geometry: geometry,
        type: type,
        tags: tags || null,
        min: [2, 1], // initial bbox values;
        max: [-1, 0]  // note that coords are usually in [0..1] range
    };
    calcBBox(feature);
    return feature;
}

function project(lonlats, tolerance) {
    var projected = [];
    for (var i = 0; i < lonlats.length; i++) {
        projected.push(projectPoint(lonlats[i]));
    }
    if (tolerance) {
        simplify(projected, tolerance);
        calcSize(projected);
    }
    return projected;
}

function projectPoint(p) {
    var sin = Math.sin(p[1] * Math.PI / 180),
        x = (p[0] / 360 + 0.5),
        y = (0.5 - 0.25 * Math.log((1 + sin) / (1 - sin)) / Math.PI);

    y = y < 0 ? 0 :
        y > 1 ? 1 : y;

    return [x, y, 0];
}

// calculate area and length of the poly
function calcSize(points) {
    var area = 0,
        dist = 0;

    for (var i = 0, a, b; i < points.length - 1; i++) {
        a = b || points[i];
        b = points[i + 1];

        area += a[0] * b[1] - b[0] * a[1];

        // use Manhattan distance instead of Euclidian one to avoid expensive square root computation
        dist += Math.abs(b[0] - a[0]) + Math.abs(b[1] - a[1]);
    }
    points.area = Math.abs(area / 2);
    points.dist = dist;
}

// calculate the feature bounding box for faster clipping later
function calcBBox(feature) {
    var geometry = feature.geometry,
        min = feature.min,
        max = feature.max;

    if (feature.type === 1) calcRingBBox(min, max, geometry);
    else for (var i = 0; i < geometry.length; i++) calcRingBBox(min, max, geometry[i]);

    return feature;
}

function calcRingBBox(min, max, points) {
    for (var i = 0, p; i < points.length; i++) {
        p = points[i];
        min[0] = Math.min(p[0], min[0]);
        max[0] = Math.max(p[0], max[0]);
        min[1] = Math.min(p[1], min[1]);
        max[1] = Math.max(p[1], max[1]);
    }
}

},{"./simplify":4}],3:[function(require,module,exports){
'use strict';

module.exports = geojsonvt;

var convert = require('./convert'),     // GeoJSON conversion and preprocessing
    transform = require('./transform'), // coordinate transformation
    clip = require('./clip'),           // stripe clipping algorithm
    wrap = require('./wrap'),           // date line processing
    createTile = require('./tile');     // final simplified tile generation


function geojsonvt(data, options) {
    return new GeoJSONVT(data, options);
}

function GeoJSONVT(data, options) {
    options = this.options = extend(Object.create(this.options), options);

    var debug = options.debug;

    if (debug) console.time('preprocess data');

    var z2 = 1 << options.maxZoom, // 2^z
        features = convert(data, options.tolerance / (z2 * options.extent));

    this.tiles = {};
    this.tileCoords = [];

    if (debug) {
        console.timeEnd('preprocess data');
        console.log('index: maxZoom: %d, maxPoints: %d', options.indexMaxZoom, options.indexMaxPoints);
        console.time('generate tiles');
        this.stats = {};
        this.total = 0;
    }

    features = wrap(features, options.buffer / options.extent, intersectX);

    // start slicing from the top tile down
    if (features.length) this.splitTile(features, 0, 0, 0);

    if (debug) {
        if (features.length) console.log('features: %d, points: %d', this.tiles[0].numFeatures, this.tiles[0].numPoints);
        console.timeEnd('generate tiles');
        console.log('tiles generated:', this.total, JSON.stringify(this.stats));
    }
}

GeoJSONVT.prototype.options = {
    maxZoom: 14,            // max zoom to preserve detail on
    indexMaxZoom: 5,        // max zoom in the tile index
    indexMaxPoints: 100000, // max number of points per tile in the tile index
    solidChildren: false,   // whether to tile solid square tiles further
    tolerance: 3,           // simplification tolerance (higher means simpler)
    extent: 4096,           // tile extent
    buffer: 64,             // tile buffer on each side
    debug: 0                // logging level (0, 1 or 2)
};

GeoJSONVT.prototype.splitTile = function (features, z, x, y, cz, cx, cy) {

    var stack = [features, z, x, y],
        options = this.options,
        debug = options.debug,
        solid = null;

    // avoid recursion by using a processing queue
    while (stack.length) {
        y = stack.pop();
        x = stack.pop();
        z = stack.pop();
        features = stack.pop();

        var z2 = 1 << z,
            id = toID(z, x, y),
            tile = this.tiles[id],
            tileTolerance = z === options.maxZoom ? 0 : options.tolerance / (z2 * options.extent);

        if (!tile) {
            if (debug > 1) console.time('creation');

            tile = this.tiles[id] = createTile(features, z2, x, y, tileTolerance, z === options.maxZoom);
            this.tileCoords.push({z: z, x: x, y: y});

            if (debug) {
                if (debug > 1) {
                    console.log('tile z%d-%d-%d (features: %d, points: %d, simplified: %d)',
                        z, x, y, tile.numFeatures, tile.numPoints, tile.numSimplified);
                    console.timeEnd('creation');
                }
                var key = 'z' + z;
                this.stats[key] = (this.stats[key] || 0) + 1;
                this.total++;
            }
        }

        // save reference to original geometry in tile so that we can drill down later if we stop now
        tile.source = features;

        // if it's the first-pass tiling
        if (!cz) {
            // stop tiling if we reached max zoom, or if the tile is too simple
            if (z === options.indexMaxZoom || tile.numPoints <= options.indexMaxPoints) continue;

        // if a drilldown to a specific tile
        } else {
            // stop tiling if we reached base zoom or our target tile zoom
            if (z === options.maxZoom || z === cz) continue;

            // stop tiling if it's not an ancestor of the target tile
            var m = 1 << (cz - z);
            if (x !== Math.floor(cx / m) || y !== Math.floor(cy / m)) continue;
        }

        // stop tiling if the tile is solid clipped square
        if (!options.solidChildren && isClippedSquare(tile, options.extent, options.buffer)) {
            if (cz) solid = z; // and remember the zoom if we're drilling down
            continue;
        }

        // if we slice further down, no need to keep source geometry
        tile.source = null;

        if (debug > 1) console.time('clipping');

        // values we'll use for clipping
        var k1 = 0.5 * options.buffer / options.extent,
            k2 = 0.5 - k1,
            k3 = 0.5 + k1,
            k4 = 1 + k1,
            tl, bl, tr, br, left, right;

        tl = bl = tr = br = null;

        left  = clip(features, z2, x - k1, x + k3, 0, intersectX, tile.min[0], tile.max[0]);
        right = clip(features, z2, x + k2, x + k4, 0, intersectX, tile.min[0], tile.max[0]);

        if (left) {
            tl = clip(left, z2, y - k1, y + k3, 1, intersectY, tile.min[1], tile.max[1]);
            bl = clip(left, z2, y + k2, y + k4, 1, intersectY, tile.min[1], tile.max[1]);
        }

        if (right) {
            tr = clip(right, z2, y - k1, y + k3, 1, intersectY, tile.min[1], tile.max[1]);
            br = clip(right, z2, y + k2, y + k4, 1, intersectY, tile.min[1], tile.max[1]);
        }

        if (debug > 1) console.timeEnd('clipping');

        if (tl) stack.push(tl, z + 1, x * 2,     y * 2);
        if (bl) stack.push(bl, z + 1, x * 2,     y * 2 + 1);
        if (tr) stack.push(tr, z + 1, x * 2 + 1, y * 2);
        if (br) stack.push(br, z + 1, x * 2 + 1, y * 2 + 1);
    }

    return solid;
};

GeoJSONVT.prototype.getTile = function (z, x, y) {
    var options = this.options,
        extent = options.extent,
        debug = options.debug;

    var z2 = 1 << z;
    x = ((x % z2) + z2) % z2; // wrap tile x coordinate

    var id = toID(z, x, y);
    if (this.tiles[id]) return transform.tile(this.tiles[id], extent);

    if (debug > 1) console.log('drilling down to z%d-%d-%d', z, x, y);

    var z0 = z,
        x0 = x,
        y0 = y,
        parent;

    while (!parent && z0 > 0) {
        z0--;
        x0 = Math.floor(x0 / 2);
        y0 = Math.floor(y0 / 2);
        parent = this.tiles[toID(z0, x0, y0)];
    }

    if (!parent || !parent.source) return null;

    // if we found a parent tile containing the original geometry, we can drill down from it
    if (debug > 1) console.log('found parent tile z%d-%d-%d', z0, x0, y0);

    // it parent tile is a solid clipped square, return it instead since it's identical
    if (isClippedSquare(parent, extent, options.buffer)) return transform.tile(parent, extent);

    if (debug > 1) console.time('drilling down');
    var solid = this.splitTile(parent.source, z0, x0, y0, z, x, y);
    if (debug > 1) console.timeEnd('drilling down');

    // one of the parent tiles was a solid clipped square
    if (solid !== null) {
        var m = 1 << (z - solid);
        id = toID(solid, Math.floor(x / m), Math.floor(y / m));
    }

    return this.tiles[id] ? transform.tile(this.tiles[id], extent) : null;
};

function toID(z, x, y) {
    return (((1 << z) * y + x) * 32) + z;
}

function intersectX(a, b, x) {
    return [x, (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0]) + a[1], 1];
}
function intersectY(a, b, y) {
    return [(y - a[1]) * (b[0] - a[0]) / (b[1] - a[1]) + a[0], y, 1];
}

function extend(dest, src) {
    for (var i in src) dest[i] = src[i];
    return dest;
}

// checks whether a tile is a whole-area fill after clipping; if it is, there's no sense slicing it further
function isClippedSquare(tile, extent, buffer) {

    var features = tile.source;
    if (features.length !== 1) return false;

    var feature = features[0];
    if (feature.type !== 3 || feature.geometry.length > 1) return false;

    var len = feature.geometry[0].length;
    if (len !== 5) return false;

    for (var i = 0; i < len; i++) {
        var p = transform.point(feature.geometry[0][i], extent, tile.z2, tile.x, tile.y);
        if ((p[0] !== -buffer && p[0] !== extent + buffer) ||
            (p[1] !== -buffer && p[1] !== extent + buffer)) return false;
    }

    return true;
}

},{"./clip":1,"./convert":2,"./tile":5,"./transform":6,"./wrap":7}],4:[function(require,module,exports){
'use strict';

module.exports = simplify;

// calculate simplification data using optimized Douglas-Peucker algorithm

function simplify(points, tolerance) {

    var sqTolerance = tolerance * tolerance,
        len = points.length,
        first = 0,
        last = len - 1,
        stack = [],
        i, maxSqDist, sqDist, index;

    // always retain the endpoints (1 is the max value)
    points[first][2] = 1;
    points[last][2] = 1;

    // avoid recursion by using a stack
    while (last) {

        maxSqDist = 0;

        for (i = first + 1; i < last; i++) {
            sqDist = getSqSegDist(points[i], points[first], points[last]);

            if (sqDist > maxSqDist) {
                index = i;
                maxSqDist = sqDist;
            }
        }

        if (maxSqDist > sqTolerance) {
            points[index][2] = maxSqDist; // save the point importance in squared pixels as a z coordinate
            stack.push(first);
            stack.push(index);
            first = index;

        } else {
            last = stack.pop();
            first = stack.pop();
        }
    }
}

// square distance from a point to a segment
function getSqSegDist(p, a, b) {

    var x = a[0], y = a[1],
        bx = b[0], by = b[1],
        px = p[0], py = p[1],
        dx = bx - x,
        dy = by - y;

    if (dx !== 0 || dy !== 0) {

        var t = ((px - x) * dx + (py - y) * dy) / (dx * dx + dy * dy);

        if (t > 1) {
            x = bx;
            y = by;

        } else if (t > 0) {
            x += dx * t;
            y += dy * t;
        }
    }

    dx = px - x;
    dy = py - y;

    return dx * dx + dy * dy;
}

},{}],5:[function(require,module,exports){
'use strict';

module.exports = createTile;

function createTile(features, z2, tx, ty, tolerance, noSimplify) {
    var tile = {
        features: [],
        numPoints: 0,
        numSimplified: 0,
        numFeatures: 0,
        source: null,
        x: tx,
        y: ty,
        z2: z2,
        transformed: false,
        min: [2, 1],
        max: [-1, 0]
    };
    for (var i = 0; i < features.length; i++) {
        tile.numFeatures++;
        addFeature(tile, features[i], tolerance, noSimplify);

        var min = features[i].min,
            max = features[i].max;

        if (min[0] < tile.min[0]) tile.min[0] = min[0];
        if (min[1] < tile.min[1]) tile.min[1] = min[1];
        if (max[0] > tile.max[0]) tile.max[0] = max[0];
        if (max[1] > tile.max[1]) tile.max[1] = max[1];
    }
    return tile;
}

function addFeature(tile, feature, tolerance, noSimplify) {

    var geom = feature.geometry,
        type = feature.type,
        simplified = [],
        sqTolerance = tolerance * tolerance,
        i, j, ring, p;

    if (type === 1) {
        for (i = 0; i < geom.length; i++) {
            simplified.push(geom[i]);
            tile.numPoints++;
            tile.numSimplified++;
        }

    } else {

        // simplify and transform projected coordinates for tile geometry
        for (i = 0; i < geom.length; i++) {
            ring = geom[i];

            // filter out tiny polylines & polygons
            if (!noSimplify && ((type === 2 && ring.dist < tolerance) ||
                                (type === 3 && ring.area < sqTolerance))) {
                tile.numPoints += ring.length;
                continue;
            }

            var simplifiedRing = [];

            for (j = 0; j < ring.length; j++) {
                p = ring[j];
                // keep points with importance > tolerance
                if (noSimplify || p[2] > sqTolerance) {
                    simplifiedRing.push(p);
                    tile.numSimplified++;
                }
                tile.numPoints++;
            }

            simplified.push(simplifiedRing);
        }
    }

    if (simplified.length) {
        tile.features.push({
            geometry: simplified,
            type: type,
            tags: feature.tags || null
        });
    }
}

},{}],6:[function(require,module,exports){
'use strict';

exports.tile = transformTile;
exports.point = transformPoint;

// Transforms the coordinates of each feature in the given tile from
// mercator-projected space into (extent x extent) tile space.
function transformTile(tile, extent) {
    if (tile.transformed) return tile;

    var z2 = tile.z2,
        tx = tile.x,
        ty = tile.y,
        i, j, k;

    for (i = 0; i < tile.features.length; i++) {
        var feature = tile.features[i],
            geom = feature.geometry,
            type = feature.type;

        if (type === 1) {
            for (j = 0; j < geom.length; j++) geom[j] = transformPoint(geom[j], extent, z2, tx, ty);

        } else {
            for (j = 0; j < geom.length; j++) {
                var ring = geom[j];
                for (k = 0; k < ring.length; k++) ring[k] = transformPoint(ring[k], extent, z2, tx, ty);
            }
        }
    }

    tile.transformed = true;

    return tile;
}

function transformPoint(p, extent, z2, tx, ty) {
    var x = Math.round(extent * (p[0] * z2 - tx)),
        y = Math.round(extent * (p[1] * z2 - ty));
    return [x, y];
}

},{}],7:[function(require,module,exports){
'use strict';

var clip = require('./clip');

module.exports = wrap;

function wrap(features, buffer, intersectX) {
    var merged = features,
        left  = clip(features, 1, -1 - buffer, buffer,     0, intersectX, -1, 2), // left world copy
        right = clip(features, 1,  1 - buffer, 2 + buffer, 0, intersectX, -1, 2); // right world copy

    if (left || right) {
        merged = clip(features, 1, -buffer, 1 + buffer, 0, intersectX, -1, 2); // center world copy

        if (left) merged = shiftFeatureCoords(left, 1).concat(merged); // merge left into center
        if (right) merged = merged.concat(shiftFeatureCoords(right, -1)); // merge right into center
    }

    return merged;
}

function shiftFeatureCoords(features, offset) {
    var newFeatures = [];

    for (var i = 0; i < features.length; i++) {
        var feature = features[i],
            type = feature.type;

        var newGeometry;

        if (type === 1) {
            newGeometry = shiftCoords(feature.geometry, offset);
        } else {
            newGeometry = [];
            for (var j = 0; j < feature.geometry.length; j++) {
                newGeometry.push(shiftCoords(feature.geometry[j], offset));
            }
        }

        newFeatures.push({
            geometry: newGeometry,
            type: type,
            tags: feature.tags,
            min: [feature.min[0] + offset, feature.min[1]],
            max: [feature.max[0] + offset, feature.max[1]]
        });
    }

    return newFeatures;
}

function shiftCoords(points, offset) {
    var newPoints = [];
    newPoints.area = points.area;
    newPoints.dist = points.dist;

    for (var i = 0; i < points.length; i++) {
        newPoints.push([points[i][0] + offset, points[i][1], points[i][2]]);
    }
    return newPoints;
}

},{"./clip":1}]},{},[3])(3)
});
//# sourceMappingURL=data:application/json;charset=utf-8;base64,eyJ2ZXJzaW9uIjozLCJzb3VyY2VzIjpbIm5vZGVfbW9kdWxlcy9icm93c2VyLXBhY2svX3ByZWx1ZGUuanMiLCJzcmMvY2xpcC5qcyIsInNyYy9jb252ZXJ0LmpzIiwic3JjL2luZGV4LmpzIiwic3JjL3NpbXBsaWZ5LmpzIiwic3JjL3RpbGUuanMiLCJzcmMvdHJhbnNmb3JtLmpzIiwic3JjL3dyYXAuanMiXSwibmFtZXMiOltdLCJtYXBwaW5ncyI6IkFBQUE7QUNBQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ3ZKQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDckpBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBOztBQ2hQQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDMUVBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7O0FDckZBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTs7QUN6Q0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQTtBQUNBO0FBQ0E7QUFDQSIsImZpbGUiOiJnZW5lcmF0ZWQuanMiLCJzb3VyY2VSb290IjoiIiwic291cmNlc0NvbnRlbnQiOlsiKGZ1bmN0aW9uIGUodCxuLHIpe2Z1bmN0aW9uIHMobyx1KXtpZighbltvXSl7aWYoIXRbb10pe3ZhciBhPXR5cGVvZiByZXF1aXJlPT1cImZ1bmN0aW9uXCImJnJlcXVpcmU7aWYoIXUmJmEpcmV0dXJuIGEobywhMCk7aWYoaSlyZXR1cm4gaShvLCEwKTt2YXIgZj1uZXcgRXJyb3IoXCJDYW5ub3QgZmluZCBtb2R1bGUgJ1wiK28rXCInXCIpO3Rocm93IGYuY29kZT1cIk1PRFVMRV9OT1RfRk9VTkRcIixmfXZhciBsPW5bb109e2V4cG9ydHM6e319O3Rbb11bMF0uY2FsbChsLmV4cG9ydHMsZnVuY3Rpb24oZSl7dmFyIG49dFtvXVsxXVtlXTtyZXR1cm4gcyhuP246ZSl9LGwsbC5leHBvcnRzLGUsdCxuLHIpfXJldHVybiBuW29dLmV4cG9ydHN9dmFyIGk9dHlwZW9mIHJlcXVpcmU9PVwiZnVuY3Rpb25cIiYmcmVxdWlyZTtmb3IodmFyIG89MDtvPHIubGVuZ3RoO28rKylzKHJbb10pO3JldHVybiBzfSkiLCIndXNlIHN0cmljdCc7XG5cbm1vZHVsZS5leHBvcnRzID0gY2xpcDtcblxuLyogY2xpcCBmZWF0dXJlcyBiZXR3ZWVuIHR3byBheGlzLXBhcmFsbGVsIGxpbmVzOlxuICogICAgIHwgICAgICAgIHxcbiAqICBfX198X19fICAgICB8ICAgICAvXG4gKiAvICAgfCAgIFxcX19fX3xfX19fL1xuICogICAgIHwgICAgICAgIHxcbiAqL1xuXG5mdW5jdGlvbiBjbGlwKGZlYXR1cmVzLCBzY2FsZSwgazEsIGsyLCBheGlzLCBpbnRlcnNlY3QsIG1pbkFsbCwgbWF4QWxsKSB7XG5cbiAgICBrMSAvPSBzY2FsZTtcbiAgICBrMiAvPSBzY2FsZTtcblxuICAgIGlmIChtaW5BbGwgPj0gazEgJiYgbWF4QWxsIDw9IGsyKSByZXR1cm4gZmVhdHVyZXM7IC8vIHRyaXZpYWwgYWNjZXB0XG4gICAgZWxzZSBpZiAobWluQWxsID4gazIgfHwgbWF4QWxsIDwgazEpIHJldHVybiBudWxsOyAvLyB0cml2aWFsIHJlamVjdFxuXG4gICAgdmFyIGNsaXBwZWQgPSBbXTtcblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZmVhdHVyZXMubGVuZ3RoOyBpKyspIHtcblxuICAgICAgICB2YXIgZmVhdHVyZSA9IGZlYXR1cmVzW2ldLFxuICAgICAgICAgICAgZ2VvbWV0cnkgPSBmZWF0dXJlLmdlb21ldHJ5LFxuICAgICAgICAgICAgdHlwZSA9IGZlYXR1cmUudHlwZSxcbiAgICAgICAgICAgIG1pbiwgbWF4O1xuXG4gICAgICAgIG1pbiA9IGZlYXR1cmUubWluW2F4aXNdO1xuICAgICAgICBtYXggPSBmZWF0dXJlLm1heFtheGlzXTtcblxuICAgICAgICBpZiAobWluID49IGsxICYmIG1heCA8PSBrMikgeyAvLyB0cml2aWFsIGFjY2VwdFxuICAgICAgICAgICAgY2xpcHBlZC5wdXNoKGZlYXR1cmUpO1xuICAgICAgICAgICAgY29udGludWU7XG4gICAgICAgIH0gZWxzZSBpZiAobWluID4gazIgfHwgbWF4IDwgazEpIGNvbnRpbnVlOyAvLyB0cml2aWFsIHJlamVjdFxuXG4gICAgICAgIHZhciBzbGljZXMgPSB0eXBlID09PSAxID9cbiAgICAgICAgICAgICAgICBjbGlwUG9pbnRzKGdlb21ldHJ5LCBrMSwgazIsIGF4aXMpIDpcbiAgICAgICAgICAgICAgICBjbGlwR2VvbWV0cnkoZ2VvbWV0cnksIGsxLCBrMiwgYXhpcywgaW50ZXJzZWN0LCB0eXBlID09PSAzKTtcblxuICAgICAgICBpZiAoc2xpY2VzLmxlbmd0aCkge1xuICAgICAgICAgICAgLy8gaWYgYSBmZWF0dXJlIGdvdCBjbGlwcGVkLCBpdCB3aWxsIGxpa2VseSBnZXQgY2xpcHBlZCBvbiB0aGUgbmV4dCB6b29tIGxldmVsIGFzIHdlbGwsXG4gICAgICAgICAgICAvLyBzbyB0aGVyZSdzIG5vIG5lZWQgdG8gcmVjYWxjdWxhdGUgYmJveGVzXG4gICAgICAgICAgICBjbGlwcGVkLnB1c2goe1xuICAgICAgICAgICAgICAgIGdlb21ldHJ5OiBzbGljZXMsXG4gICAgICAgICAgICAgICAgdHlwZTogdHlwZSxcbiAgICAgICAgICAgICAgICB0YWdzOiBmZWF0dXJlc1tpXS50YWdzIHx8IG51bGwsXG4gICAgICAgICAgICAgICAgbWluOiBmZWF0dXJlLm1pbixcbiAgICAgICAgICAgICAgICBtYXg6IGZlYXR1cmUubWF4XG4gICAgICAgICAgICB9KTtcbiAgICAgICAgfVxuICAgIH1cblxuICAgIHJldHVybiBjbGlwcGVkLmxlbmd0aCA/IGNsaXBwZWQgOiBudWxsO1xufVxuXG5mdW5jdGlvbiBjbGlwUG9pbnRzKGdlb21ldHJ5LCBrMSwgazIsIGF4aXMpIHtcbiAgICB2YXIgc2xpY2UgPSBbXTtcblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZ2VvbWV0cnkubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgdmFyIGEgPSBnZW9tZXRyeVtpXSxcbiAgICAgICAgICAgIGFrID0gYVtheGlzXTtcblxuICAgICAgICBpZiAoYWsgPj0gazEgJiYgYWsgPD0gazIpIHNsaWNlLnB1c2goYSk7XG4gICAgfVxuICAgIHJldHVybiBzbGljZTtcbn1cblxuZnVuY3Rpb24gY2xpcEdlb21ldHJ5KGdlb21ldHJ5LCBrMSwgazIsIGF4aXMsIGludGVyc2VjdCwgY2xvc2VkKSB7XG5cbiAgICB2YXIgc2xpY2VzID0gW107XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGdlb21ldHJ5Lmxlbmd0aDsgaSsrKSB7XG5cbiAgICAgICAgdmFyIGFrID0gMCxcbiAgICAgICAgICAgIGJrID0gMCxcbiAgICAgICAgICAgIGIgPSBudWxsLFxuICAgICAgICAgICAgcG9pbnRzID0gZ2VvbWV0cnlbaV0sXG4gICAgICAgICAgICBhcmVhID0gcG9pbnRzLmFyZWEsXG4gICAgICAgICAgICBkaXN0ID0gcG9pbnRzLmRpc3QsXG4gICAgICAgICAgICBsZW4gPSBwb2ludHMubGVuZ3RoLFxuICAgICAgICAgICAgYSwgaiwgbGFzdDtcblxuICAgICAgICB2YXIgc2xpY2UgPSBbXTtcblxuICAgICAgICBmb3IgKGogPSAwOyBqIDwgbGVuIC0gMTsgaisrKSB7XG4gICAgICAgICAgICBhID0gYiB8fCBwb2ludHNbal07XG4gICAgICAgICAgICBiID0gcG9pbnRzW2ogKyAxXTtcbiAgICAgICAgICAgIGFrID0gYmsgfHwgYVtheGlzXTtcbiAgICAgICAgICAgIGJrID0gYltheGlzXTtcblxuICAgICAgICAgICAgaWYgKGFrIDwgazEpIHtcblxuICAgICAgICAgICAgICAgIGlmICgoYmsgPiBrMikpIHsgLy8gLS0tfC0tLS0tfC0tPlxuICAgICAgICAgICAgICAgICAgICBzbGljZS5wdXNoKGludGVyc2VjdChhLCBiLCBrMSksIGludGVyc2VjdChhLCBiLCBrMikpO1xuICAgICAgICAgICAgICAgICAgICBpZiAoIWNsb3NlZCkgc2xpY2UgPSBuZXdTbGljZShzbGljZXMsIHNsaWNlLCBhcmVhLCBkaXN0KTtcblxuICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAoYmsgPj0gazEpIHNsaWNlLnB1c2goaW50ZXJzZWN0KGEsIGIsIGsxKSk7IC8vIC0tLXwtLT4gIHxcblxuICAgICAgICAgICAgfSBlbHNlIGlmIChhayA+IGsyKSB7XG5cbiAgICAgICAgICAgICAgICBpZiAoKGJrIDwgazEpKSB7IC8vIDwtLXwtLS0tLXwtLS1cbiAgICAgICAgICAgICAgICAgICAgc2xpY2UucHVzaChpbnRlcnNlY3QoYSwgYiwgazIpLCBpbnRlcnNlY3QoYSwgYiwgazEpKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKCFjbG9zZWQpIHNsaWNlID0gbmV3U2xpY2Uoc2xpY2VzLCBzbGljZSwgYXJlYSwgZGlzdCk7XG5cbiAgICAgICAgICAgICAgICB9IGVsc2UgaWYgKGJrIDw9IGsyKSBzbGljZS5wdXNoKGludGVyc2VjdChhLCBiLCBrMikpOyAvLyB8ICA8LS18LS0tXG5cbiAgICAgICAgICAgIH0gZWxzZSB7XG5cbiAgICAgICAgICAgICAgICBzbGljZS5wdXNoKGEpO1xuXG4gICAgICAgICAgICAgICAgaWYgKGJrIDwgazEpIHsgLy8gPC0tfC0tLSAgfFxuICAgICAgICAgICAgICAgICAgICBzbGljZS5wdXNoKGludGVyc2VjdChhLCBiLCBrMSkpO1xuICAgICAgICAgICAgICAgICAgICBpZiAoIWNsb3NlZCkgc2xpY2UgPSBuZXdTbGljZShzbGljZXMsIHNsaWNlLCBhcmVhLCBkaXN0KTtcblxuICAgICAgICAgICAgICAgIH0gZWxzZSBpZiAoYmsgPiBrMikgeyAvLyB8ICAtLS18LS0+XG4gICAgICAgICAgICAgICAgICAgIHNsaWNlLnB1c2goaW50ZXJzZWN0KGEsIGIsIGsyKSk7XG4gICAgICAgICAgICAgICAgICAgIGlmICghY2xvc2VkKSBzbGljZSA9IG5ld1NsaWNlKHNsaWNlcywgc2xpY2UsIGFyZWEsIGRpc3QpO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgICAvLyB8IC0tPiB8XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cblxuICAgICAgICAvLyBhZGQgdGhlIGxhc3QgcG9pbnRcbiAgICAgICAgYSA9IHBvaW50c1tsZW4gLSAxXTtcbiAgICAgICAgYWsgPSBhW2F4aXNdO1xuICAgICAgICBpZiAoYWsgPj0gazEgJiYgYWsgPD0gazIpIHNsaWNlLnB1c2goYSk7XG5cbiAgICAgICAgLy8gY2xvc2UgdGhlIHBvbHlnb24gaWYgaXRzIGVuZHBvaW50cyBhcmUgbm90IHRoZSBzYW1lIGFmdGVyIGNsaXBwaW5nXG5cbiAgICAgICAgbGFzdCA9IHNsaWNlW3NsaWNlLmxlbmd0aCAtIDFdO1xuICAgICAgICBpZiAoY2xvc2VkICYmIGxhc3QgJiYgKHNsaWNlWzBdWzBdICE9PSBsYXN0WzBdIHx8IHNsaWNlWzBdWzFdICE9PSBsYXN0WzFdKSkgc2xpY2UucHVzaChzbGljZVswXSk7XG5cbiAgICAgICAgLy8gYWRkIHRoZSBmaW5hbCBzbGljZVxuICAgICAgICBuZXdTbGljZShzbGljZXMsIHNsaWNlLCBhcmVhLCBkaXN0KTtcbiAgICB9XG5cbiAgICByZXR1cm4gc2xpY2VzO1xufVxuXG5mdW5jdGlvbiBuZXdTbGljZShzbGljZXMsIHNsaWNlLCBhcmVhLCBkaXN0KSB7XG4gICAgaWYgKHNsaWNlLmxlbmd0aCkge1xuICAgICAgICAvLyB3ZSBkb24ndCByZWNhbGN1bGF0ZSB0aGUgYXJlYS9sZW5ndGggb2YgdGhlIHVuY2xpcHBlZCBnZW9tZXRyeSBiZWNhdXNlIHRoZSBjYXNlIHdoZXJlIGl0IGdvZXNcbiAgICAgICAgLy8gYmVsb3cgdGhlIHZpc2liaWxpdHkgdGhyZXNob2xkIGFzIGEgcmVzdWx0IG9mIGNsaXBwaW5nIGlzIHJhcmUsIHNvIHdlIGF2b2lkIGRvaW5nIHVubmVjZXNzYXJ5IHdvcmtcbiAgICAgICAgc2xpY2UuYXJlYSA9IGFyZWE7XG4gICAgICAgIHNsaWNlLmRpc3QgPSBkaXN0O1xuXG4gICAgICAgIHNsaWNlcy5wdXNoKHNsaWNlKTtcbiAgICB9XG4gICAgcmV0dXJuIFtdO1xufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG5tb2R1bGUuZXhwb3J0cyA9IGNvbnZlcnQ7XG5cbnZhciBzaW1wbGlmeSA9IHJlcXVpcmUoJy4vc2ltcGxpZnknKTtcblxuLy8gY29udmVydHMgR2VvSlNPTiBmZWF0dXJlIGludG8gYW4gaW50ZXJtZWRpYXRlIHByb2plY3RlZCBKU09OIHZlY3RvciBmb3JtYXQgd2l0aCBzaW1wbGlmaWNhdGlvbiBkYXRhXG5cbmZ1bmN0aW9uIGNvbnZlcnQoZGF0YSwgdG9sZXJhbmNlKSB7XG4gICAgdmFyIGZlYXR1cmVzID0gW107XG5cbiAgICBpZiAoZGF0YS50eXBlID09PSAnRmVhdHVyZUNvbGxlY3Rpb24nKSB7XG4gICAgICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZGF0YS5mZWF0dXJlcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgY29udmVydEZlYXR1cmUoZmVhdHVyZXMsIGRhdGEuZmVhdHVyZXNbaV0sIHRvbGVyYW5jZSk7XG4gICAgICAgIH1cbiAgICB9IGVsc2UgaWYgKGRhdGEudHlwZSA9PT0gJ0ZlYXR1cmUnKSB7XG4gICAgICAgIGNvbnZlcnRGZWF0dXJlKGZlYXR1cmVzLCBkYXRhLCB0b2xlcmFuY2UpO1xuXG4gICAgfSBlbHNlIHtcbiAgICAgICAgLy8gc2luZ2xlIGdlb21ldHJ5IG9yIGEgZ2VvbWV0cnkgY29sbGVjdGlvblxuICAgICAgICBjb252ZXJ0RmVhdHVyZShmZWF0dXJlcywge2dlb21ldHJ5OiBkYXRhfSwgdG9sZXJhbmNlKTtcbiAgICB9XG4gICAgcmV0dXJuIGZlYXR1cmVzO1xufVxuXG5mdW5jdGlvbiBjb252ZXJ0RmVhdHVyZShmZWF0dXJlcywgZmVhdHVyZSwgdG9sZXJhbmNlKSB7XG4gICAgaWYgKGZlYXR1cmUuZ2VvbWV0cnkgPT09IG51bGwpIHtcbiAgICAgICAgLy8gaWdub3JlIGZlYXR1cmVzIHdpdGggbnVsbCBnZW9tZXRyeVxuICAgICAgICByZXR1cm47XG4gICAgfVxuXG4gICAgdmFyIGdlb20gPSBmZWF0dXJlLmdlb21ldHJ5LFxuICAgICAgICB0eXBlID0gZ2VvbS50eXBlLFxuICAgICAgICBjb29yZHMgPSBnZW9tLmNvb3JkaW5hdGVzLFxuICAgICAgICB0YWdzID0gZmVhdHVyZS5wcm9wZXJ0aWVzLFxuICAgICAgICBpLCBqLCByaW5ncztcblxuICAgIGlmICh0eXBlID09PSAnUG9pbnQnKSB7XG4gICAgICAgIGZlYXR1cmVzLnB1c2goY3JlYXRlKHRhZ3MsIDEsIFtwcm9qZWN0UG9pbnQoY29vcmRzKV0pKTtcblxuICAgIH0gZWxzZSBpZiAodHlwZSA9PT0gJ011bHRpUG9pbnQnKSB7XG4gICAgICAgIGZlYXR1cmVzLnB1c2goY3JlYXRlKHRhZ3MsIDEsIHByb2plY3QoY29vcmRzKSkpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnTGluZVN0cmluZycpIHtcbiAgICAgICAgZmVhdHVyZXMucHVzaChjcmVhdGUodGFncywgMiwgW3Byb2plY3QoY29vcmRzLCB0b2xlcmFuY2UpXSkpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnTXVsdGlMaW5lU3RyaW5nJyB8fCB0eXBlID09PSAnUG9seWdvbicpIHtcbiAgICAgICAgcmluZ3MgPSBbXTtcbiAgICAgICAgZm9yIChpID0gMDsgaSA8IGNvb3Jkcy5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgcmluZ3MucHVzaChwcm9qZWN0KGNvb3Jkc1tpXSwgdG9sZXJhbmNlKSk7XG4gICAgICAgIH1cbiAgICAgICAgZmVhdHVyZXMucHVzaChjcmVhdGUodGFncywgdHlwZSA9PT0gJ1BvbHlnb24nID8gMyA6IDIsIHJpbmdzKSk7XG5cbiAgICB9IGVsc2UgaWYgKHR5cGUgPT09ICdNdWx0aVBvbHlnb24nKSB7XG4gICAgICAgIHJpbmdzID0gW107XG4gICAgICAgIGZvciAoaSA9IDA7IGkgPCBjb29yZHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgICAgIGZvciAoaiA9IDA7IGogPCBjb29yZHNbaV0ubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICByaW5ncy5wdXNoKHByb2plY3QoY29vcmRzW2ldW2pdLCB0b2xlcmFuY2UpKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuICAgICAgICBmZWF0dXJlcy5wdXNoKGNyZWF0ZSh0YWdzLCAzLCByaW5ncykpO1xuXG4gICAgfSBlbHNlIGlmICh0eXBlID09PSAnR2VvbWV0cnlDb2xsZWN0aW9uJykge1xuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZ2VvbS5nZW9tZXRyaWVzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgICAgICBjb252ZXJ0RmVhdHVyZShmZWF0dXJlcywge1xuICAgICAgICAgICAgICAgIGdlb21ldHJ5OiBnZW9tLmdlb21ldHJpZXNbaV0sXG4gICAgICAgICAgICAgICAgcHJvcGVydGllczogdGFnc1xuICAgICAgICAgICAgfSwgdG9sZXJhbmNlKTtcbiAgICAgICAgfVxuXG4gICAgfSBlbHNlIHtcbiAgICAgICAgdGhyb3cgbmV3IEVycm9yKCdJbnB1dCBkYXRhIGlzIG5vdCBhIHZhbGlkIEdlb0pTT04gb2JqZWN0LicpO1xuICAgIH1cbn1cblxuZnVuY3Rpb24gY3JlYXRlKHRhZ3MsIHR5cGUsIGdlb21ldHJ5KSB7XG4gICAgdmFyIGZlYXR1cmUgPSB7XG4gICAgICAgIGdlb21ldHJ5OiBnZW9tZXRyeSxcbiAgICAgICAgdHlwZTogdHlwZSxcbiAgICAgICAgdGFnczogdGFncyB8fCBudWxsLFxuICAgICAgICBtaW46IFsyLCAxXSwgLy8gaW5pdGlhbCBiYm94IHZhbHVlcztcbiAgICAgICAgbWF4OiBbLTEsIDBdICAvLyBub3RlIHRoYXQgY29vcmRzIGFyZSB1c3VhbGx5IGluIFswLi4xXSByYW5nZVxuICAgIH07XG4gICAgY2FsY0JCb3goZmVhdHVyZSk7XG4gICAgcmV0dXJuIGZlYXR1cmU7XG59XG5cbmZ1bmN0aW9uIHByb2plY3QobG9ubGF0cywgdG9sZXJhbmNlKSB7XG4gICAgdmFyIHByb2plY3RlZCA9IFtdO1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgbG9ubGF0cy5sZW5ndGg7IGkrKykge1xuICAgICAgICBwcm9qZWN0ZWQucHVzaChwcm9qZWN0UG9pbnQobG9ubGF0c1tpXSkpO1xuICAgIH1cbiAgICBpZiAodG9sZXJhbmNlKSB7XG4gICAgICAgIHNpbXBsaWZ5KHByb2plY3RlZCwgdG9sZXJhbmNlKTtcbiAgICAgICAgY2FsY1NpemUocHJvamVjdGVkKTtcbiAgICB9XG4gICAgcmV0dXJuIHByb2plY3RlZDtcbn1cblxuZnVuY3Rpb24gcHJvamVjdFBvaW50KHApIHtcbiAgICB2YXIgc2luID0gTWF0aC5zaW4ocFsxXSAqIE1hdGguUEkgLyAxODApLFxuICAgICAgICB4ID0gKHBbMF0gLyAzNjAgKyAwLjUpLFxuICAgICAgICB5ID0gKDAuNSAtIDAuMjUgKiBNYXRoLmxvZygoMSArIHNpbikgLyAoMSAtIHNpbikpIC8gTWF0aC5QSSk7XG5cbiAgICB5ID0geSA8IDAgPyAwIDpcbiAgICAgICAgeSA+IDEgPyAxIDogeTtcblxuICAgIHJldHVybiBbeCwgeSwgMF07XG59XG5cbi8vIGNhbGN1bGF0ZSBhcmVhIGFuZCBsZW5ndGggb2YgdGhlIHBvbHlcbmZ1bmN0aW9uIGNhbGNTaXplKHBvaW50cykge1xuICAgIHZhciBhcmVhID0gMCxcbiAgICAgICAgZGlzdCA9IDA7XG5cbiAgICBmb3IgKHZhciBpID0gMCwgYSwgYjsgaSA8IHBvaW50cy5sZW5ndGggLSAxOyBpKyspIHtcbiAgICAgICAgYSA9IGIgfHwgcG9pbnRzW2ldO1xuICAgICAgICBiID0gcG9pbnRzW2kgKyAxXTtcblxuICAgICAgICBhcmVhICs9IGFbMF0gKiBiWzFdIC0gYlswXSAqIGFbMV07XG5cbiAgICAgICAgLy8gdXNlIE1hbmhhdHRhbiBkaXN0YW5jZSBpbnN0ZWFkIG9mIEV1Y2xpZGlhbiBvbmUgdG8gYXZvaWQgZXhwZW5zaXZlIHNxdWFyZSByb290IGNvbXB1dGF0aW9uXG4gICAgICAgIGRpc3QgKz0gTWF0aC5hYnMoYlswXSAtIGFbMF0pICsgTWF0aC5hYnMoYlsxXSAtIGFbMV0pO1xuICAgIH1cbiAgICBwb2ludHMuYXJlYSA9IE1hdGguYWJzKGFyZWEgLyAyKTtcbiAgICBwb2ludHMuZGlzdCA9IGRpc3Q7XG59XG5cbi8vIGNhbGN1bGF0ZSB0aGUgZmVhdHVyZSBib3VuZGluZyBib3ggZm9yIGZhc3RlciBjbGlwcGluZyBsYXRlclxuZnVuY3Rpb24gY2FsY0JCb3goZmVhdHVyZSkge1xuICAgIHZhciBnZW9tZXRyeSA9IGZlYXR1cmUuZ2VvbWV0cnksXG4gICAgICAgIG1pbiA9IGZlYXR1cmUubWluLFxuICAgICAgICBtYXggPSBmZWF0dXJlLm1heDtcblxuICAgIGlmIChmZWF0dXJlLnR5cGUgPT09IDEpIGNhbGNSaW5nQkJveChtaW4sIG1heCwgZ2VvbWV0cnkpO1xuICAgIGVsc2UgZm9yICh2YXIgaSA9IDA7IGkgPCBnZW9tZXRyeS5sZW5ndGg7IGkrKykgY2FsY1JpbmdCQm94KG1pbiwgbWF4LCBnZW9tZXRyeVtpXSk7XG5cbiAgICByZXR1cm4gZmVhdHVyZTtcbn1cblxuZnVuY3Rpb24gY2FsY1JpbmdCQm94KG1pbiwgbWF4LCBwb2ludHMpIHtcbiAgICBmb3IgKHZhciBpID0gMCwgcDsgaSA8IHBvaW50cy5sZW5ndGg7IGkrKykge1xuICAgICAgICBwID0gcG9pbnRzW2ldO1xuICAgICAgICBtaW5bMF0gPSBNYXRoLm1pbihwWzBdLCBtaW5bMF0pO1xuICAgICAgICBtYXhbMF0gPSBNYXRoLm1heChwWzBdLCBtYXhbMF0pO1xuICAgICAgICBtaW5bMV0gPSBNYXRoLm1pbihwWzFdLCBtaW5bMV0pO1xuICAgICAgICBtYXhbMV0gPSBNYXRoLm1heChwWzFdLCBtYXhbMV0pO1xuICAgIH1cbn1cbiIsIid1c2Ugc3RyaWN0JztcblxubW9kdWxlLmV4cG9ydHMgPSBnZW9qc29udnQ7XG5cbnZhciBjb252ZXJ0ID0gcmVxdWlyZSgnLi9jb252ZXJ0JyksICAgICAvLyBHZW9KU09OIGNvbnZlcnNpb24gYW5kIHByZXByb2Nlc3NpbmdcbiAgICB0cmFuc2Zvcm0gPSByZXF1aXJlKCcuL3RyYW5zZm9ybScpLCAvLyBjb29yZGluYXRlIHRyYW5zZm9ybWF0aW9uXG4gICAgY2xpcCA9IHJlcXVpcmUoJy4vY2xpcCcpLCAgICAgICAgICAgLy8gc3RyaXBlIGNsaXBwaW5nIGFsZ29yaXRobVxuICAgIHdyYXAgPSByZXF1aXJlKCcuL3dyYXAnKSwgICAgICAgICAgIC8vIGRhdGUgbGluZSBwcm9jZXNzaW5nXG4gICAgY3JlYXRlVGlsZSA9IHJlcXVpcmUoJy4vdGlsZScpOyAgICAgLy8gZmluYWwgc2ltcGxpZmllZCB0aWxlIGdlbmVyYXRpb25cblxuXG5mdW5jdGlvbiBnZW9qc29udnQoZGF0YSwgb3B0aW9ucykge1xuICAgIHJldHVybiBuZXcgR2VvSlNPTlZUKGRhdGEsIG9wdGlvbnMpO1xufVxuXG5mdW5jdGlvbiBHZW9KU09OVlQoZGF0YSwgb3B0aW9ucykge1xuICAgIG9wdGlvbnMgPSB0aGlzLm9wdGlvbnMgPSBleHRlbmQoT2JqZWN0LmNyZWF0ZSh0aGlzLm9wdGlvbnMpLCBvcHRpb25zKTtcblxuICAgIHZhciBkZWJ1ZyA9IG9wdGlvbnMuZGVidWc7XG5cbiAgICBpZiAoZGVidWcpIGNvbnNvbGUudGltZSgncHJlcHJvY2VzcyBkYXRhJyk7XG5cbiAgICB2YXIgejIgPSAxIDw8IG9wdGlvbnMubWF4Wm9vbSwgLy8gMl56XG4gICAgICAgIGZlYXR1cmVzID0gY29udmVydChkYXRhLCBvcHRpb25zLnRvbGVyYW5jZSAvICh6MiAqIG9wdGlvbnMuZXh0ZW50KSk7XG5cbiAgICB0aGlzLnRpbGVzID0ge307XG4gICAgdGhpcy50aWxlQ29vcmRzID0gW107XG5cbiAgICBpZiAoZGVidWcpIHtcbiAgICAgICAgY29uc29sZS50aW1lRW5kKCdwcmVwcm9jZXNzIGRhdGEnKTtcbiAgICAgICAgY29uc29sZS5sb2coJ2luZGV4OiBtYXhab29tOiAlZCwgbWF4UG9pbnRzOiAlZCcsIG9wdGlvbnMuaW5kZXhNYXhab29tLCBvcHRpb25zLmluZGV4TWF4UG9pbnRzKTtcbiAgICAgICAgY29uc29sZS50aW1lKCdnZW5lcmF0ZSB0aWxlcycpO1xuICAgICAgICB0aGlzLnN0YXRzID0ge307XG4gICAgICAgIHRoaXMudG90YWwgPSAwO1xuICAgIH1cblxuICAgIGZlYXR1cmVzID0gd3JhcChmZWF0dXJlcywgb3B0aW9ucy5idWZmZXIgLyBvcHRpb25zLmV4dGVudCwgaW50ZXJzZWN0WCk7XG5cbiAgICAvLyBzdGFydCBzbGljaW5nIGZyb20gdGhlIHRvcCB0aWxlIGRvd25cbiAgICBpZiAoZmVhdHVyZXMubGVuZ3RoKSB0aGlzLnNwbGl0VGlsZShmZWF0dXJlcywgMCwgMCwgMCk7XG5cbiAgICBpZiAoZGVidWcpIHtcbiAgICAgICAgaWYgKGZlYXR1cmVzLmxlbmd0aCkgY29uc29sZS5sb2coJ2ZlYXR1cmVzOiAlZCwgcG9pbnRzOiAlZCcsIHRoaXMudGlsZXNbMF0ubnVtRmVhdHVyZXMsIHRoaXMudGlsZXNbMF0ubnVtUG9pbnRzKTtcbiAgICAgICAgY29uc29sZS50aW1lRW5kKCdnZW5lcmF0ZSB0aWxlcycpO1xuICAgICAgICBjb25zb2xlLmxvZygndGlsZXMgZ2VuZXJhdGVkOicsIHRoaXMudG90YWwsIEpTT04uc3RyaW5naWZ5KHRoaXMuc3RhdHMpKTtcbiAgICB9XG59XG5cbkdlb0pTT05WVC5wcm90b3R5cGUub3B0aW9ucyA9IHtcbiAgICBtYXhab29tOiAxNCwgICAgICAgICAgICAvLyBtYXggem9vbSB0byBwcmVzZXJ2ZSBkZXRhaWwgb25cbiAgICBpbmRleE1heFpvb206IDUsICAgICAgICAvLyBtYXggem9vbSBpbiB0aGUgdGlsZSBpbmRleFxuICAgIGluZGV4TWF4UG9pbnRzOiAxMDAwMDAsIC8vIG1heCBudW1iZXIgb2YgcG9pbnRzIHBlciB0aWxlIGluIHRoZSB0aWxlIGluZGV4XG4gICAgc29saWRDaGlsZHJlbjogZmFsc2UsICAgLy8gd2hldGhlciB0byB0aWxlIHNvbGlkIHNxdWFyZSB0aWxlcyBmdXJ0aGVyXG4gICAgdG9sZXJhbmNlOiAzLCAgICAgICAgICAgLy8gc2ltcGxpZmljYXRpb24gdG9sZXJhbmNlIChoaWdoZXIgbWVhbnMgc2ltcGxlcilcbiAgICBleHRlbnQ6IDQwOTYsICAgICAgICAgICAvLyB0aWxlIGV4dGVudFxuICAgIGJ1ZmZlcjogNjQsICAgICAgICAgICAgIC8vIHRpbGUgYnVmZmVyIG9uIGVhY2ggc2lkZVxuICAgIGRlYnVnOiAwICAgICAgICAgICAgICAgIC8vIGxvZ2dpbmcgbGV2ZWwgKDAsIDEgb3IgMilcbn07XG5cbkdlb0pTT05WVC5wcm90b3R5cGUuc3BsaXRUaWxlID0gZnVuY3Rpb24gKGZlYXR1cmVzLCB6LCB4LCB5LCBjeiwgY3gsIGN5KSB7XG5cbiAgICB2YXIgc3RhY2sgPSBbZmVhdHVyZXMsIHosIHgsIHldLFxuICAgICAgICBvcHRpb25zID0gdGhpcy5vcHRpb25zLFxuICAgICAgICBkZWJ1ZyA9IG9wdGlvbnMuZGVidWcsXG4gICAgICAgIHNvbGlkID0gbnVsbDtcblxuICAgIC8vIGF2b2lkIHJlY3Vyc2lvbiBieSB1c2luZyBhIHByb2Nlc3NpbmcgcXVldWVcbiAgICB3aGlsZSAoc3RhY2subGVuZ3RoKSB7XG4gICAgICAgIHkgPSBzdGFjay5wb3AoKTtcbiAgICAgICAgeCA9IHN0YWNrLnBvcCgpO1xuICAgICAgICB6ID0gc3RhY2sucG9wKCk7XG4gICAgICAgIGZlYXR1cmVzID0gc3RhY2sucG9wKCk7XG5cbiAgICAgICAgdmFyIHoyID0gMSA8PCB6LFxuICAgICAgICAgICAgaWQgPSB0b0lEKHosIHgsIHkpLFxuICAgICAgICAgICAgdGlsZSA9IHRoaXMudGlsZXNbaWRdLFxuICAgICAgICAgICAgdGlsZVRvbGVyYW5jZSA9IHogPT09IG9wdGlvbnMubWF4Wm9vbSA/IDAgOiBvcHRpb25zLnRvbGVyYW5jZSAvICh6MiAqIG9wdGlvbnMuZXh0ZW50KTtcblxuICAgICAgICBpZiAoIXRpbGUpIHtcbiAgICAgICAgICAgIGlmIChkZWJ1ZyA+IDEpIGNvbnNvbGUudGltZSgnY3JlYXRpb24nKTtcblxuICAgICAgICAgICAgdGlsZSA9IHRoaXMudGlsZXNbaWRdID0gY3JlYXRlVGlsZShmZWF0dXJlcywgejIsIHgsIHksIHRpbGVUb2xlcmFuY2UsIHogPT09IG9wdGlvbnMubWF4Wm9vbSk7XG4gICAgICAgICAgICB0aGlzLnRpbGVDb29yZHMucHVzaCh7ejogeiwgeDogeCwgeTogeX0pO1xuXG4gICAgICAgICAgICBpZiAoZGVidWcpIHtcbiAgICAgICAgICAgICAgICBpZiAoZGVidWcgPiAxKSB7XG4gICAgICAgICAgICAgICAgICAgIGNvbnNvbGUubG9nKCd0aWxlIHolZC0lZC0lZCAoZmVhdHVyZXM6ICVkLCBwb2ludHM6ICVkLCBzaW1wbGlmaWVkOiAlZCknLFxuICAgICAgICAgICAgICAgICAgICAgICAgeiwgeCwgeSwgdGlsZS5udW1GZWF0dXJlcywgdGlsZS5udW1Qb2ludHMsIHRpbGUubnVtU2ltcGxpZmllZCk7XG4gICAgICAgICAgICAgICAgICAgIGNvbnNvbGUudGltZUVuZCgnY3JlYXRpb24nKTtcbiAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgdmFyIGtleSA9ICd6JyArIHo7XG4gICAgICAgICAgICAgICAgdGhpcy5zdGF0c1trZXldID0gKHRoaXMuc3RhdHNba2V5XSB8fCAwKSArIDE7XG4gICAgICAgICAgICAgICAgdGhpcy50b3RhbCsrO1xuICAgICAgICAgICAgfVxuICAgICAgICB9XG5cbiAgICAgICAgLy8gc2F2ZSByZWZlcmVuY2UgdG8gb3JpZ2luYWwgZ2VvbWV0cnkgaW4gdGlsZSBzbyB0aGF0IHdlIGNhbiBkcmlsbCBkb3duIGxhdGVyIGlmIHdlIHN0b3Agbm93XG4gICAgICAgIHRpbGUuc291cmNlID0gZmVhdHVyZXM7XG5cbiAgICAgICAgLy8gaWYgaXQncyB0aGUgZmlyc3QtcGFzcyB0aWxpbmdcbiAgICAgICAgaWYgKCFjeikge1xuICAgICAgICAgICAgLy8gc3RvcCB0aWxpbmcgaWYgd2UgcmVhY2hlZCBtYXggem9vbSwgb3IgaWYgdGhlIHRpbGUgaXMgdG9vIHNpbXBsZVxuICAgICAgICAgICAgaWYgKHogPT09IG9wdGlvbnMuaW5kZXhNYXhab29tIHx8IHRpbGUubnVtUG9pbnRzIDw9IG9wdGlvbnMuaW5kZXhNYXhQb2ludHMpIGNvbnRpbnVlO1xuXG4gICAgICAgIC8vIGlmIGEgZHJpbGxkb3duIHRvIGEgc3BlY2lmaWMgdGlsZVxuICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgLy8gc3RvcCB0aWxpbmcgaWYgd2UgcmVhY2hlZCBiYXNlIHpvb20gb3Igb3VyIHRhcmdldCB0aWxlIHpvb21cbiAgICAgICAgICAgIGlmICh6ID09PSBvcHRpb25zLm1heFpvb20gfHwgeiA9PT0gY3opIGNvbnRpbnVlO1xuXG4gICAgICAgICAgICAvLyBzdG9wIHRpbGluZyBpZiBpdCdzIG5vdCBhbiBhbmNlc3RvciBvZiB0aGUgdGFyZ2V0IHRpbGVcbiAgICAgICAgICAgIHZhciBtID0gMSA8PCAoY3ogLSB6KTtcbiAgICAgICAgICAgIGlmICh4ICE9PSBNYXRoLmZsb29yKGN4IC8gbSkgfHwgeSAhPT0gTWF0aC5mbG9vcihjeSAvIG0pKSBjb250aW51ZTtcbiAgICAgICAgfVxuXG4gICAgICAgIC8vIHN0b3AgdGlsaW5nIGlmIHRoZSB0aWxlIGlzIHNvbGlkIGNsaXBwZWQgc3F1YXJlXG4gICAgICAgIGlmICghb3B0aW9ucy5zb2xpZENoaWxkcmVuICYmIGlzQ2xpcHBlZFNxdWFyZSh0aWxlLCBvcHRpb25zLmV4dGVudCwgb3B0aW9ucy5idWZmZXIpKSB7XG4gICAgICAgICAgICBpZiAoY3opIHNvbGlkID0gejsgLy8gYW5kIHJlbWVtYmVyIHRoZSB6b29tIGlmIHdlJ3JlIGRyaWxsaW5nIGRvd25cbiAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICB9XG5cbiAgICAgICAgLy8gaWYgd2Ugc2xpY2UgZnVydGhlciBkb3duLCBubyBuZWVkIHRvIGtlZXAgc291cmNlIGdlb21ldHJ5XG4gICAgICAgIHRpbGUuc291cmNlID0gbnVsbDtcblxuICAgICAgICBpZiAoZGVidWcgPiAxKSBjb25zb2xlLnRpbWUoJ2NsaXBwaW5nJyk7XG5cbiAgICAgICAgLy8gdmFsdWVzIHdlJ2xsIHVzZSBmb3IgY2xpcHBpbmdcbiAgICAgICAgdmFyIGsxID0gMC41ICogb3B0aW9ucy5idWZmZXIgLyBvcHRpb25zLmV4dGVudCxcbiAgICAgICAgICAgIGsyID0gMC41IC0gazEsXG4gICAgICAgICAgICBrMyA9IDAuNSArIGsxLFxuICAgICAgICAgICAgazQgPSAxICsgazEsXG4gICAgICAgICAgICB0bCwgYmwsIHRyLCBiciwgbGVmdCwgcmlnaHQ7XG5cbiAgICAgICAgdGwgPSBibCA9IHRyID0gYnIgPSBudWxsO1xuXG4gICAgICAgIGxlZnQgID0gY2xpcChmZWF0dXJlcywgejIsIHggLSBrMSwgeCArIGszLCAwLCBpbnRlcnNlY3RYLCB0aWxlLm1pblswXSwgdGlsZS5tYXhbMF0pO1xuICAgICAgICByaWdodCA9IGNsaXAoZmVhdHVyZXMsIHoyLCB4ICsgazIsIHggKyBrNCwgMCwgaW50ZXJzZWN0WCwgdGlsZS5taW5bMF0sIHRpbGUubWF4WzBdKTtcblxuICAgICAgICBpZiAobGVmdCkge1xuICAgICAgICAgICAgdGwgPSBjbGlwKGxlZnQsIHoyLCB5IC0gazEsIHkgKyBrMywgMSwgaW50ZXJzZWN0WSwgdGlsZS5taW5bMV0sIHRpbGUubWF4WzFdKTtcbiAgICAgICAgICAgIGJsID0gY2xpcChsZWZ0LCB6MiwgeSArIGsyLCB5ICsgazQsIDEsIGludGVyc2VjdFksIHRpbGUubWluWzFdLCB0aWxlLm1heFsxXSk7XG4gICAgICAgIH1cblxuICAgICAgICBpZiAocmlnaHQpIHtcbiAgICAgICAgICAgIHRyID0gY2xpcChyaWdodCwgejIsIHkgLSBrMSwgeSArIGszLCAxLCBpbnRlcnNlY3RZLCB0aWxlLm1pblsxXSwgdGlsZS5tYXhbMV0pO1xuICAgICAgICAgICAgYnIgPSBjbGlwKHJpZ2h0LCB6MiwgeSArIGsyLCB5ICsgazQsIDEsIGludGVyc2VjdFksIHRpbGUubWluWzFdLCB0aWxlLm1heFsxXSk7XG4gICAgICAgIH1cblxuICAgICAgICBpZiAoZGVidWcgPiAxKSBjb25zb2xlLnRpbWVFbmQoJ2NsaXBwaW5nJyk7XG5cbiAgICAgICAgaWYgKHRsKSBzdGFjay5wdXNoKHRsLCB6ICsgMSwgeCAqIDIsICAgICB5ICogMik7XG4gICAgICAgIGlmIChibCkgc3RhY2sucHVzaChibCwgeiArIDEsIHggKiAyLCAgICAgeSAqIDIgKyAxKTtcbiAgICAgICAgaWYgKHRyKSBzdGFjay5wdXNoKHRyLCB6ICsgMSwgeCAqIDIgKyAxLCB5ICogMik7XG4gICAgICAgIGlmIChicikgc3RhY2sucHVzaChiciwgeiArIDEsIHggKiAyICsgMSwgeSAqIDIgKyAxKTtcbiAgICB9XG5cbiAgICByZXR1cm4gc29saWQ7XG59O1xuXG5HZW9KU09OVlQucHJvdG90eXBlLmdldFRpbGUgPSBmdW5jdGlvbiAoeiwgeCwgeSkge1xuICAgIHZhciBvcHRpb25zID0gdGhpcy5vcHRpb25zLFxuICAgICAgICBleHRlbnQgPSBvcHRpb25zLmV4dGVudCxcbiAgICAgICAgZGVidWcgPSBvcHRpb25zLmRlYnVnO1xuXG4gICAgdmFyIHoyID0gMSA8PCB6O1xuICAgIHggPSAoKHggJSB6MikgKyB6MikgJSB6MjsgLy8gd3JhcCB0aWxlIHggY29vcmRpbmF0ZVxuXG4gICAgdmFyIGlkID0gdG9JRCh6LCB4LCB5KTtcbiAgICBpZiAodGhpcy50aWxlc1tpZF0pIHJldHVybiB0cmFuc2Zvcm0udGlsZSh0aGlzLnRpbGVzW2lkXSwgZXh0ZW50KTtcblxuICAgIGlmIChkZWJ1ZyA+IDEpIGNvbnNvbGUubG9nKCdkcmlsbGluZyBkb3duIHRvIHolZC0lZC0lZCcsIHosIHgsIHkpO1xuXG4gICAgdmFyIHowID0geixcbiAgICAgICAgeDAgPSB4LFxuICAgICAgICB5MCA9IHksXG4gICAgICAgIHBhcmVudDtcblxuICAgIHdoaWxlICghcGFyZW50ICYmIHowID4gMCkge1xuICAgICAgICB6MC0tO1xuICAgICAgICB4MCA9IE1hdGguZmxvb3IoeDAgLyAyKTtcbiAgICAgICAgeTAgPSBNYXRoLmZsb29yKHkwIC8gMik7XG4gICAgICAgIHBhcmVudCA9IHRoaXMudGlsZXNbdG9JRCh6MCwgeDAsIHkwKV07XG4gICAgfVxuXG4gICAgaWYgKCFwYXJlbnQgfHwgIXBhcmVudC5zb3VyY2UpIHJldHVybiBudWxsO1xuXG4gICAgLy8gaWYgd2UgZm91bmQgYSBwYXJlbnQgdGlsZSBjb250YWluaW5nIHRoZSBvcmlnaW5hbCBnZW9tZXRyeSwgd2UgY2FuIGRyaWxsIGRvd24gZnJvbSBpdFxuICAgIGlmIChkZWJ1ZyA+IDEpIGNvbnNvbGUubG9nKCdmb3VuZCBwYXJlbnQgdGlsZSB6JWQtJWQtJWQnLCB6MCwgeDAsIHkwKTtcblxuICAgIC8vIGl0IHBhcmVudCB0aWxlIGlzIGEgc29saWQgY2xpcHBlZCBzcXVhcmUsIHJldHVybiBpdCBpbnN0ZWFkIHNpbmNlIGl0J3MgaWRlbnRpY2FsXG4gICAgaWYgKGlzQ2xpcHBlZFNxdWFyZShwYXJlbnQsIGV4dGVudCwgb3B0aW9ucy5idWZmZXIpKSByZXR1cm4gdHJhbnNmb3JtLnRpbGUocGFyZW50LCBleHRlbnQpO1xuXG4gICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS50aW1lKCdkcmlsbGluZyBkb3duJyk7XG4gICAgdmFyIHNvbGlkID0gdGhpcy5zcGxpdFRpbGUocGFyZW50LnNvdXJjZSwgejAsIHgwLCB5MCwgeiwgeCwgeSk7XG4gICAgaWYgKGRlYnVnID4gMSkgY29uc29sZS50aW1lRW5kKCdkcmlsbGluZyBkb3duJyk7XG5cbiAgICAvLyBvbmUgb2YgdGhlIHBhcmVudCB0aWxlcyB3YXMgYSBzb2xpZCBjbGlwcGVkIHNxdWFyZVxuICAgIGlmIChzb2xpZCAhPT0gbnVsbCkge1xuICAgICAgICB2YXIgbSA9IDEgPDwgKHogLSBzb2xpZCk7XG4gICAgICAgIGlkID0gdG9JRChzb2xpZCwgTWF0aC5mbG9vcih4IC8gbSksIE1hdGguZmxvb3IoeSAvIG0pKTtcbiAgICB9XG5cbiAgICByZXR1cm4gdGhpcy50aWxlc1tpZF0gPyB0cmFuc2Zvcm0udGlsZSh0aGlzLnRpbGVzW2lkXSwgZXh0ZW50KSA6IG51bGw7XG59O1xuXG5mdW5jdGlvbiB0b0lEKHosIHgsIHkpIHtcbiAgICByZXR1cm4gKCgoMSA8PCB6KSAqIHkgKyB4KSAqIDMyKSArIHo7XG59XG5cbmZ1bmN0aW9uIGludGVyc2VjdFgoYSwgYiwgeCkge1xuICAgIHJldHVybiBbeCwgKHggLSBhWzBdKSAqIChiWzFdIC0gYVsxXSkgLyAoYlswXSAtIGFbMF0pICsgYVsxXSwgMV07XG59XG5mdW5jdGlvbiBpbnRlcnNlY3RZKGEsIGIsIHkpIHtcbiAgICByZXR1cm4gWyh5IC0gYVsxXSkgKiAoYlswXSAtIGFbMF0pIC8gKGJbMV0gLSBhWzFdKSArIGFbMF0sIHksIDFdO1xufVxuXG5mdW5jdGlvbiBleHRlbmQoZGVzdCwgc3JjKSB7XG4gICAgZm9yICh2YXIgaSBpbiBzcmMpIGRlc3RbaV0gPSBzcmNbaV07XG4gICAgcmV0dXJuIGRlc3Q7XG59XG5cbi8vIGNoZWNrcyB3aGV0aGVyIGEgdGlsZSBpcyBhIHdob2xlLWFyZWEgZmlsbCBhZnRlciBjbGlwcGluZzsgaWYgaXQgaXMsIHRoZXJlJ3Mgbm8gc2Vuc2Ugc2xpY2luZyBpdCBmdXJ0aGVyXG5mdW5jdGlvbiBpc0NsaXBwZWRTcXVhcmUodGlsZSwgZXh0ZW50LCBidWZmZXIpIHtcblxuICAgIHZhciBmZWF0dXJlcyA9IHRpbGUuc291cmNlO1xuICAgIGlmIChmZWF0dXJlcy5sZW5ndGggIT09IDEpIHJldHVybiBmYWxzZTtcblxuICAgIHZhciBmZWF0dXJlID0gZmVhdHVyZXNbMF07XG4gICAgaWYgKGZlYXR1cmUudHlwZSAhPT0gMyB8fCBmZWF0dXJlLmdlb21ldHJ5Lmxlbmd0aCA+IDEpIHJldHVybiBmYWxzZTtcblxuICAgIHZhciBsZW4gPSBmZWF0dXJlLmdlb21ldHJ5WzBdLmxlbmd0aDtcbiAgICBpZiAobGVuICE9PSA1KSByZXR1cm4gZmFsc2U7XG5cbiAgICBmb3IgKHZhciBpID0gMDsgaSA8IGxlbjsgaSsrKSB7XG4gICAgICAgIHZhciBwID0gdHJhbnNmb3JtLnBvaW50KGZlYXR1cmUuZ2VvbWV0cnlbMF1baV0sIGV4dGVudCwgdGlsZS56MiwgdGlsZS54LCB0aWxlLnkpO1xuICAgICAgICBpZiAoKHBbMF0gIT09IC1idWZmZXIgJiYgcFswXSAhPT0gZXh0ZW50ICsgYnVmZmVyKSB8fFxuICAgICAgICAgICAgKHBbMV0gIT09IC1idWZmZXIgJiYgcFsxXSAhPT0gZXh0ZW50ICsgYnVmZmVyKSkgcmV0dXJuIGZhbHNlO1xuICAgIH1cblxuICAgIHJldHVybiB0cnVlO1xufVxuIiwiJ3VzZSBzdHJpY3QnO1xuXG5tb2R1bGUuZXhwb3J0cyA9IHNpbXBsaWZ5O1xuXG4vLyBjYWxjdWxhdGUgc2ltcGxpZmljYXRpb24gZGF0YSB1c2luZyBvcHRpbWl6ZWQgRG91Z2xhcy1QZXVja2VyIGFsZ29yaXRobVxuXG5mdW5jdGlvbiBzaW1wbGlmeShwb2ludHMsIHRvbGVyYW5jZSkge1xuXG4gICAgdmFyIHNxVG9sZXJhbmNlID0gdG9sZXJhbmNlICogdG9sZXJhbmNlLFxuICAgICAgICBsZW4gPSBwb2ludHMubGVuZ3RoLFxuICAgICAgICBmaXJzdCA9IDAsXG4gICAgICAgIGxhc3QgPSBsZW4gLSAxLFxuICAgICAgICBzdGFjayA9IFtdLFxuICAgICAgICBpLCBtYXhTcURpc3QsIHNxRGlzdCwgaW5kZXg7XG5cbiAgICAvLyBhbHdheXMgcmV0YWluIHRoZSBlbmRwb2ludHMgKDEgaXMgdGhlIG1heCB2YWx1ZSlcbiAgICBwb2ludHNbZmlyc3RdWzJdID0gMTtcbiAgICBwb2ludHNbbGFzdF1bMl0gPSAxO1xuXG4gICAgLy8gYXZvaWQgcmVjdXJzaW9uIGJ5IHVzaW5nIGEgc3RhY2tcbiAgICB3aGlsZSAobGFzdCkge1xuXG4gICAgICAgIG1heFNxRGlzdCA9IDA7XG5cbiAgICAgICAgZm9yIChpID0gZmlyc3QgKyAxOyBpIDwgbGFzdDsgaSsrKSB7XG4gICAgICAgICAgICBzcURpc3QgPSBnZXRTcVNlZ0Rpc3QocG9pbnRzW2ldLCBwb2ludHNbZmlyc3RdLCBwb2ludHNbbGFzdF0pO1xuXG4gICAgICAgICAgICBpZiAoc3FEaXN0ID4gbWF4U3FEaXN0KSB7XG4gICAgICAgICAgICAgICAgaW5kZXggPSBpO1xuICAgICAgICAgICAgICAgIG1heFNxRGlzdCA9IHNxRGlzdDtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIGlmIChtYXhTcURpc3QgPiBzcVRvbGVyYW5jZSkge1xuICAgICAgICAgICAgcG9pbnRzW2luZGV4XVsyXSA9IG1heFNxRGlzdDsgLy8gc2F2ZSB0aGUgcG9pbnQgaW1wb3J0YW5jZSBpbiBzcXVhcmVkIHBpeGVscyBhcyBhIHogY29vcmRpbmF0ZVxuICAgICAgICAgICAgc3RhY2sucHVzaChmaXJzdCk7XG4gICAgICAgICAgICBzdGFjay5wdXNoKGluZGV4KTtcbiAgICAgICAgICAgIGZpcnN0ID0gaW5kZXg7XG5cbiAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGxhc3QgPSBzdGFjay5wb3AoKTtcbiAgICAgICAgICAgIGZpcnN0ID0gc3RhY2sucG9wKCk7XG4gICAgICAgIH1cbiAgICB9XG59XG5cbi8vIHNxdWFyZSBkaXN0YW5jZSBmcm9tIGEgcG9pbnQgdG8gYSBzZWdtZW50XG5mdW5jdGlvbiBnZXRTcVNlZ0Rpc3QocCwgYSwgYikge1xuXG4gICAgdmFyIHggPSBhWzBdLCB5ID0gYVsxXSxcbiAgICAgICAgYnggPSBiWzBdLCBieSA9IGJbMV0sXG4gICAgICAgIHB4ID0gcFswXSwgcHkgPSBwWzFdLFxuICAgICAgICBkeCA9IGJ4IC0geCxcbiAgICAgICAgZHkgPSBieSAtIHk7XG5cbiAgICBpZiAoZHggIT09IDAgfHwgZHkgIT09IDApIHtcblxuICAgICAgICB2YXIgdCA9ICgocHggLSB4KSAqIGR4ICsgKHB5IC0geSkgKiBkeSkgLyAoZHggKiBkeCArIGR5ICogZHkpO1xuXG4gICAgICAgIGlmICh0ID4gMSkge1xuICAgICAgICAgICAgeCA9IGJ4O1xuICAgICAgICAgICAgeSA9IGJ5O1xuXG4gICAgICAgIH0gZWxzZSBpZiAodCA+IDApIHtcbiAgICAgICAgICAgIHggKz0gZHggKiB0O1xuICAgICAgICAgICAgeSArPSBkeSAqIHQ7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBkeCA9IHB4IC0geDtcbiAgICBkeSA9IHB5IC0geTtcblxuICAgIHJldHVybiBkeCAqIGR4ICsgZHkgKiBkeTtcbn1cbiIsIid1c2Ugc3RyaWN0JztcblxubW9kdWxlLmV4cG9ydHMgPSBjcmVhdGVUaWxlO1xuXG5mdW5jdGlvbiBjcmVhdGVUaWxlKGZlYXR1cmVzLCB6MiwgdHgsIHR5LCB0b2xlcmFuY2UsIG5vU2ltcGxpZnkpIHtcbiAgICB2YXIgdGlsZSA9IHtcbiAgICAgICAgZmVhdHVyZXM6IFtdLFxuICAgICAgICBudW1Qb2ludHM6IDAsXG4gICAgICAgIG51bVNpbXBsaWZpZWQ6IDAsXG4gICAgICAgIG51bUZlYXR1cmVzOiAwLFxuICAgICAgICBzb3VyY2U6IG51bGwsXG4gICAgICAgIHg6IHR4LFxuICAgICAgICB5OiB0eSxcbiAgICAgICAgejI6IHoyLFxuICAgICAgICB0cmFuc2Zvcm1lZDogZmFsc2UsXG4gICAgICAgIG1pbjogWzIsIDFdLFxuICAgICAgICBtYXg6IFstMSwgMF1cbiAgICB9O1xuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZmVhdHVyZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgdGlsZS5udW1GZWF0dXJlcysrO1xuICAgICAgICBhZGRGZWF0dXJlKHRpbGUsIGZlYXR1cmVzW2ldLCB0b2xlcmFuY2UsIG5vU2ltcGxpZnkpO1xuXG4gICAgICAgIHZhciBtaW4gPSBmZWF0dXJlc1tpXS5taW4sXG4gICAgICAgICAgICBtYXggPSBmZWF0dXJlc1tpXS5tYXg7XG5cbiAgICAgICAgaWYgKG1pblswXSA8IHRpbGUubWluWzBdKSB0aWxlLm1pblswXSA9IG1pblswXTtcbiAgICAgICAgaWYgKG1pblsxXSA8IHRpbGUubWluWzFdKSB0aWxlLm1pblsxXSA9IG1pblsxXTtcbiAgICAgICAgaWYgKG1heFswXSA+IHRpbGUubWF4WzBdKSB0aWxlLm1heFswXSA9IG1heFswXTtcbiAgICAgICAgaWYgKG1heFsxXSA+IHRpbGUubWF4WzFdKSB0aWxlLm1heFsxXSA9IG1heFsxXTtcbiAgICB9XG4gICAgcmV0dXJuIHRpbGU7XG59XG5cbmZ1bmN0aW9uIGFkZEZlYXR1cmUodGlsZSwgZmVhdHVyZSwgdG9sZXJhbmNlLCBub1NpbXBsaWZ5KSB7XG5cbiAgICB2YXIgZ2VvbSA9IGZlYXR1cmUuZ2VvbWV0cnksXG4gICAgICAgIHR5cGUgPSBmZWF0dXJlLnR5cGUsXG4gICAgICAgIHNpbXBsaWZpZWQgPSBbXSxcbiAgICAgICAgc3FUb2xlcmFuY2UgPSB0b2xlcmFuY2UgKiB0b2xlcmFuY2UsXG4gICAgICAgIGksIGosIHJpbmcsIHA7XG5cbiAgICBpZiAodHlwZSA9PT0gMSkge1xuICAgICAgICBmb3IgKGkgPSAwOyBpIDwgZ2VvbS5sZW5ndGg7IGkrKykge1xuICAgICAgICAgICAgc2ltcGxpZmllZC5wdXNoKGdlb21baV0pO1xuICAgICAgICAgICAgdGlsZS5udW1Qb2ludHMrKztcbiAgICAgICAgICAgIHRpbGUubnVtU2ltcGxpZmllZCsrO1xuICAgICAgICB9XG5cbiAgICB9IGVsc2Uge1xuXG4gICAgICAgIC8vIHNpbXBsaWZ5IGFuZCB0cmFuc2Zvcm0gcHJvamVjdGVkIGNvb3JkaW5hdGVzIGZvciB0aWxlIGdlb21ldHJ5XG4gICAgICAgIGZvciAoaSA9IDA7IGkgPCBnZW9tLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgICAgICByaW5nID0gZ2VvbVtpXTtcblxuICAgICAgICAgICAgLy8gZmlsdGVyIG91dCB0aW55IHBvbHlsaW5lcyAmIHBvbHlnb25zXG4gICAgICAgICAgICBpZiAoIW5vU2ltcGxpZnkgJiYgKCh0eXBlID09PSAyICYmIHJpbmcuZGlzdCA8IHRvbGVyYW5jZSkgfHxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgKHR5cGUgPT09IDMgJiYgcmluZy5hcmVhIDwgc3FUb2xlcmFuY2UpKSkge1xuICAgICAgICAgICAgICAgIHRpbGUubnVtUG9pbnRzICs9IHJpbmcubGVuZ3RoO1xuICAgICAgICAgICAgICAgIGNvbnRpbnVlO1xuICAgICAgICAgICAgfVxuXG4gICAgICAgICAgICB2YXIgc2ltcGxpZmllZFJpbmcgPSBbXTtcblxuICAgICAgICAgICAgZm9yIChqID0gMDsgaiA8IHJpbmcubGVuZ3RoOyBqKyspIHtcbiAgICAgICAgICAgICAgICBwID0gcmluZ1tqXTtcbiAgICAgICAgICAgICAgICAvLyBrZWVwIHBvaW50cyB3aXRoIGltcG9ydGFuY2UgPiB0b2xlcmFuY2VcbiAgICAgICAgICAgICAgICBpZiAobm9TaW1wbGlmeSB8fCBwWzJdID4gc3FUb2xlcmFuY2UpIHtcbiAgICAgICAgICAgICAgICAgICAgc2ltcGxpZmllZFJpbmcucHVzaChwKTtcbiAgICAgICAgICAgICAgICAgICAgdGlsZS5udW1TaW1wbGlmaWVkKys7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgIHRpbGUubnVtUG9pbnRzKys7XG4gICAgICAgICAgICB9XG5cbiAgICAgICAgICAgIHNpbXBsaWZpZWQucHVzaChzaW1wbGlmaWVkUmluZyk7XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICBpZiAoc2ltcGxpZmllZC5sZW5ndGgpIHtcbiAgICAgICAgdGlsZS5mZWF0dXJlcy5wdXNoKHtcbiAgICAgICAgICAgIGdlb21ldHJ5OiBzaW1wbGlmaWVkLFxuICAgICAgICAgICAgdHlwZTogdHlwZSxcbiAgICAgICAgICAgIHRhZ3M6IGZlYXR1cmUudGFncyB8fCBudWxsXG4gICAgICAgIH0pO1xuICAgIH1cbn1cbiIsIid1c2Ugc3RyaWN0JztcblxuZXhwb3J0cy50aWxlID0gdHJhbnNmb3JtVGlsZTtcbmV4cG9ydHMucG9pbnQgPSB0cmFuc2Zvcm1Qb2ludDtcblxuLy8gVHJhbnNmb3JtcyB0aGUgY29vcmRpbmF0ZXMgb2YgZWFjaCBmZWF0dXJlIGluIHRoZSBnaXZlbiB0aWxlIGZyb21cbi8vIG1lcmNhdG9yLXByb2plY3RlZCBzcGFjZSBpbnRvIChleHRlbnQgeCBleHRlbnQpIHRpbGUgc3BhY2UuXG5mdW5jdGlvbiB0cmFuc2Zvcm1UaWxlKHRpbGUsIGV4dGVudCkge1xuICAgIGlmICh0aWxlLnRyYW5zZm9ybWVkKSByZXR1cm4gdGlsZTtcblxuICAgIHZhciB6MiA9IHRpbGUuejIsXG4gICAgICAgIHR4ID0gdGlsZS54LFxuICAgICAgICB0eSA9IHRpbGUueSxcbiAgICAgICAgaSwgaiwgaztcblxuICAgIGZvciAoaSA9IDA7IGkgPCB0aWxlLmZlYXR1cmVzLmxlbmd0aDsgaSsrKSB7XG4gICAgICAgIHZhciBmZWF0dXJlID0gdGlsZS5mZWF0dXJlc1tpXSxcbiAgICAgICAgICAgIGdlb20gPSBmZWF0dXJlLmdlb21ldHJ5LFxuICAgICAgICAgICAgdHlwZSA9IGZlYXR1cmUudHlwZTtcblxuICAgICAgICBpZiAodHlwZSA9PT0gMSkge1xuICAgICAgICAgICAgZm9yIChqID0gMDsgaiA8IGdlb20ubGVuZ3RoOyBqKyspIGdlb21bal0gPSB0cmFuc2Zvcm1Qb2ludChnZW9tW2pdLCBleHRlbnQsIHoyLCB0eCwgdHkpO1xuXG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBmb3IgKGogPSAwOyBqIDwgZ2VvbS5sZW5ndGg7IGorKykge1xuICAgICAgICAgICAgICAgIHZhciByaW5nID0gZ2VvbVtqXTtcbiAgICAgICAgICAgICAgICBmb3IgKGsgPSAwOyBrIDwgcmluZy5sZW5ndGg7IGsrKykgcmluZ1trXSA9IHRyYW5zZm9ybVBvaW50KHJpbmdba10sIGV4dGVudCwgejIsIHR4LCB0eSk7XG4gICAgICAgICAgICB9XG4gICAgICAgIH1cbiAgICB9XG5cbiAgICB0aWxlLnRyYW5zZm9ybWVkID0gdHJ1ZTtcblxuICAgIHJldHVybiB0aWxlO1xufVxuXG5mdW5jdGlvbiB0cmFuc2Zvcm1Qb2ludChwLCBleHRlbnQsIHoyLCB0eCwgdHkpIHtcbiAgICB2YXIgeCA9IE1hdGgucm91bmQoZXh0ZW50ICogKHBbMF0gKiB6MiAtIHR4KSksXG4gICAgICAgIHkgPSBNYXRoLnJvdW5kKGV4dGVudCAqIChwWzFdICogejIgLSB0eSkpO1xuICAgIHJldHVybiBbeCwgeV07XG59XG4iLCIndXNlIHN0cmljdCc7XG5cbnZhciBjbGlwID0gcmVxdWlyZSgnLi9jbGlwJyk7XG5cbm1vZHVsZS5leHBvcnRzID0gd3JhcDtcblxuZnVuY3Rpb24gd3JhcChmZWF0dXJlcywgYnVmZmVyLCBpbnRlcnNlY3RYKSB7XG4gICAgdmFyIG1lcmdlZCA9IGZlYXR1cmVzLFxuICAgICAgICBsZWZ0ICA9IGNsaXAoZmVhdHVyZXMsIDEsIC0xIC0gYnVmZmVyLCBidWZmZXIsICAgICAwLCBpbnRlcnNlY3RYLCAtMSwgMiksIC8vIGxlZnQgd29ybGQgY29weVxuICAgICAgICByaWdodCA9IGNsaXAoZmVhdHVyZXMsIDEsICAxIC0gYnVmZmVyLCAyICsgYnVmZmVyLCAwLCBpbnRlcnNlY3RYLCAtMSwgMik7IC8vIHJpZ2h0IHdvcmxkIGNvcHlcblxuICAgIGlmIChsZWZ0IHx8IHJpZ2h0KSB7XG4gICAgICAgIG1lcmdlZCA9IGNsaXAoZmVhdHVyZXMsIDEsIC1idWZmZXIsIDEgKyBidWZmZXIsIDAsIGludGVyc2VjdFgsIC0xLCAyKTsgLy8gY2VudGVyIHdvcmxkIGNvcHlcblxuICAgICAgICBpZiAobGVmdCkgbWVyZ2VkID0gc2hpZnRGZWF0dXJlQ29vcmRzKGxlZnQsIDEpLmNvbmNhdChtZXJnZWQpOyAvLyBtZXJnZSBsZWZ0IGludG8gY2VudGVyXG4gICAgICAgIGlmIChyaWdodCkgbWVyZ2VkID0gbWVyZ2VkLmNvbmNhdChzaGlmdEZlYXR1cmVDb29yZHMocmlnaHQsIC0xKSk7IC8vIG1lcmdlIHJpZ2h0IGludG8gY2VudGVyXG4gICAgfVxuXG4gICAgcmV0dXJuIG1lcmdlZDtcbn1cblxuZnVuY3Rpb24gc2hpZnRGZWF0dXJlQ29vcmRzKGZlYXR1cmVzLCBvZmZzZXQpIHtcbiAgICB2YXIgbmV3RmVhdHVyZXMgPSBbXTtcblxuICAgIGZvciAodmFyIGkgPSAwOyBpIDwgZmVhdHVyZXMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgdmFyIGZlYXR1cmUgPSBmZWF0dXJlc1tpXSxcbiAgICAgICAgICAgIHR5cGUgPSBmZWF0dXJlLnR5cGU7XG5cbiAgICAgICAgdmFyIG5ld0dlb21ldHJ5O1xuXG4gICAgICAgIGlmICh0eXBlID09PSAxKSB7XG4gICAgICAgICAgICBuZXdHZW9tZXRyeSA9IHNoaWZ0Q29vcmRzKGZlYXR1cmUuZ2VvbWV0cnksIG9mZnNldCk7XG4gICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICBuZXdHZW9tZXRyeSA9IFtdO1xuICAgICAgICAgICAgZm9yICh2YXIgaiA9IDA7IGogPCBmZWF0dXJlLmdlb21ldHJ5Lmxlbmd0aDsgaisrKSB7XG4gICAgICAgICAgICAgICAgbmV3R2VvbWV0cnkucHVzaChzaGlmdENvb3JkcyhmZWF0dXJlLmdlb21ldHJ5W2pdLCBvZmZzZXQpKTtcbiAgICAgICAgICAgIH1cbiAgICAgICAgfVxuXG4gICAgICAgIG5ld0ZlYXR1cmVzLnB1c2goe1xuICAgICAgICAgICAgZ2VvbWV0cnk6IG5ld0dlb21ldHJ5LFxuICAgICAgICAgICAgdHlwZTogdHlwZSxcbiAgICAgICAgICAgIHRhZ3M6IGZlYXR1cmUudGFncyxcbiAgICAgICAgICAgIG1pbjogW2ZlYXR1cmUubWluWzBdICsgb2Zmc2V0LCBmZWF0dXJlLm1pblsxXV0sXG4gICAgICAgICAgICBtYXg6IFtmZWF0dXJlLm1heFswXSArIG9mZnNldCwgZmVhdHVyZS5tYXhbMV1dXG4gICAgICAgIH0pO1xuICAgIH1cblxuICAgIHJldHVybiBuZXdGZWF0dXJlcztcbn1cblxuZnVuY3Rpb24gc2hpZnRDb29yZHMocG9pbnRzLCBvZmZzZXQpIHtcbiAgICB2YXIgbmV3UG9pbnRzID0gW107XG4gICAgbmV3UG9pbnRzLmFyZWEgPSBwb2ludHMuYXJlYTtcbiAgICBuZXdQb2ludHMuZGlzdCA9IHBvaW50cy5kaXN0O1xuXG4gICAgZm9yICh2YXIgaSA9IDA7IGkgPCBwb2ludHMubGVuZ3RoOyBpKyspIHtcbiAgICAgICAgbmV3UG9pbnRzLnB1c2goW3BvaW50c1tpXVswXSArIG9mZnNldCwgcG9pbnRzW2ldWzFdLCBwb2ludHNbaV1bMl1dKTtcbiAgICB9XG4gICAgcmV0dXJuIG5ld1BvaW50cztcbn1cbiJdfQ==
