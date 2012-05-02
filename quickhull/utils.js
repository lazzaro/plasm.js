function getRandomPoints3D(numPoint, xMax, yMax, zMax) {
    var points = new Array();
    var phase = Math.random() * Math.PI * 2;
    for (var i = 0; i < numPoint/2; i++) {
        var r =  Math.random()*xMax/4;
        var theta = Math.random() * 1.5 * Math.PI + phase;
        points.push( [ xMax /4 + r * Math.cos(theta), yMax/2 + 2 * r * Math.sin(theta),  zMax/5 + 2 * r * Math.sin(theta)] );
    }
    var phase = Math.random() * Math.PI * 2;
    for (var i = 0; i < numPoint/2; i++) {
        var r =  Math.random()*xMax/4;
        var theta = Math.random() * 1.5 * Math.PI + phase;
        points.push( [ xMax /4 * 3 +  r * Math.cos(theta), yMax/2 +  r * Math.sin(theta),  zMax /5 * 3 +  r * Math.cos(theta) ] );
    }
    return points;
}

function getRandomPoints2D(numPoint, xMax, yMax) {
    var points = new Array();
    var phase = Math.random() * Math.PI * 2;
    for (var i = 0; i < numPoint/2; i++) {
        var r =  Math.random()*xMax/4;
        var theta = Math.random() * 1.5 * Math.PI + phase;
        points.push( [ xMax /4 + r * Math.cos(theta), yMax/2 + 2 * r * Math.sin(theta)] );
    }
    var phase = Math.random() * Math.PI * 2;
    for (var i = 0; i < numPoint/2; i++) {
        var r =  Math.random()*xMax/4;
        var theta = Math.random() * 1.5 * Math.PI + phase;
        points.push( [ xMax /4 * 3 +  r * Math.cos(theta), yMax/2 +  r * Math.sin(theta)] );
    }
    return points;
}

function random2D(numPoints, maxX, maxY) {
	var points = new Array();

	for ( var i = 0; i < numPoints; i++) {
		
		points.push([Math.floor(Math.random() * maxX) * Math.pow(-1, Math.round(Math.random())), 
		             Math.floor(Math.random() * maxY) * Math.pow(-1, Math.round(Math.random()))]);
	}
	
	return points;
}

function random3D(numPoints, maxX, maxY, maxZ) {
	var points = new Array();

	for ( var i = 0; i < numPoints; i++) {
		
		points.push([Math.random() * maxX * Math.pow(-1, Math.round(Math.random())), 
		             Math.random() * maxY * Math.pow(-1, Math.round(Math.random())),
		             Math.random() * maxZ * Math.pow(-1, Math.round(Math.random()))]);
	}
	
	return points;
}

function qhPlotPoints(pts) {
    ctx = document.getElementById('qh_demo').getContext('2d');
    ctx.clearRect(0,0,400,400);
    ctx.fillStyle = 'rgb(0,0,0)';
    for (var idx in pts) {
        var pt = pts[idx];
        ctx.fillRect(pt[0],pt[1],2,2);
    }
}

function qhPlotPoint(pt) {
    ctx = document.getElementById('qh_demo').getContext('2d');
    ctx.fillStyle = 'rgb(255,0,0)';
    
    ctx.fillRect(pt[0],pt[1],2,2);
    
}

function plotBaseLine(baseLine,color) {
    var ctx = document.getElementById('qh_demo').getContext('2d');
    var pt1 = baseLine[0];
    var pt2 = baseLine[1];
    ctx.save();
    ctx.strokeStyle = color;
    ctx.beginPath();
    ctx.moveTo(pt1[0],pt1[1]);
    ctx.lineTo(pt2[0],pt2[1]);
    ctx.closePath();
    ctx.stroke();
    ctx.restore();
} 