//var points = getRandomPoints3D(3, 100, 100, 100);
var points = [[148, 41], [150, 299], [36, 211], [338, 237], [316, 178], [338, 207], [250, 250]];//getRandomPoints2D(50,399,399);
//[[71, 169], [44, 190], [56, 154], [360, 216], [333, 240], [340, 154]];//[[10, 10], [390, 10], [200, 350], [200, 340], [250, 340]];//
console.log(points);

var convexHull = quickhull.quickhull(points);
var vertices;

for ( var i = 0; i < convexHull.length; i++) {
	vertices = convexHull[i].vertices;
	plotBaseLine(vertices,'rgb(255,0,0)');
}