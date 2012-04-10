/*!
i * Quickhull algoritm in javascript translated from 
 * http://www.qhull.org
 * 
 */
!(function (exports) {

	/**
	 * Variables.
	 */

	var cos = Math.cos;
	var sin = Math.sin;
	var round = Math.round;
	var min = Math.min;
	var abs = Math.abs;
	var pi = Math.PI;
	var random = Math.random;
	var floor = Math.floor;

	/**
	 * Library namespace.
	 */

	var quickhull = exports.quickhull = {};

	/**
	 * Library version.
	 */

	quickhull.version = '0.0.0';

	/**
	 * utils namespace
	 * @api private
	 */

	quickhull._utils = {};



	/**
	 * verify if a point is already inside the set
	 * 
	 * @param {Array|Float32Array} set - set of points
	 * @param {Number|Float32} point - point to verify
	 * @return {Boolean} - if the point is already in
	 * 
	 */

	var isIn = quickhull._utils.isIn = function(set, point) {
		var isIn = false, equal;
		var tempPoint;

		for ( var i = 0; i < set.length && !isIn; i++) {
			tempPoint = set[i];
			equal = true;
			for ( var j = 0; j < tempPoint.length && equal; j++) {
				// se almeno una coordiata è differente non sono lo stesso punto
				if (tempPoint[j] !== point[j]) {
					equal = false;
				}
			}

			if (!equal) {
				isIn = true;
			}
		}

		return isIn;
	};

	/**
	 * Add a point to a set only if not already present
	 * 
	 * @param {Array|Float32Array} set - set of points
	 * @param {Array|Float32} point - point to add
	 * 
	 */

	var addNoDup = quickhull._utils.addNoDup = function (set, point){

		if (!isIn(set, point)) {
			set.push(point);
		}

		return set;
	};

	/**
	 * Compute the determinant of a simplex
	 * 
	 * @param {Array|Float32} point - point apex
	 * @param {Array|Float32Array} simplex - base points
	 * @param {Number} dim - dimension
	 */

	var detSimplex = quickhull._utils.detSimplex = function(point, simplex, dim) {
		
		var det2 = function(a1, a2, b1, b2) {
			return (a1 * b2) - (a2 * b1);
		};
		
		var det3 = function(a1, a2, a3, b1, b2, b3, c1, c2, c3) {
			return (a1) * det2(b2,b3,c2,c3) - (b1) * det2(a2,a3,c2,c3) + (c1)*det2(a2,a3,b2,b3);
		};
		
		var j = 0, det;
		var sign;
		var rows = new Array(), rowsGauss = new Array(), coordPointSimpl, row = new Array();
		
		for ( var i = 0; i < simplex.length && !(j == dim); i++) {
			coordPointSimpl = simplex[i];
			for ( var k = 0; k < dim; k++) {
				row.push(coordPointSimpl[k] - point[k]);
			}
			
			rows.push(row);
			row = new Array();
			j = rows.length;
		}
		
		
		// calcolo del determinante
		//TODO implementare gestione nearzero
		if (dim == 2) {
			det = det2(rows[0][0], rows[0][1], rows[1][0], rows[1][1]); 
		} else if (dim == 3) {
			det = det3(rows[0][0], rows[0][1], rows[0][2], rows[1][0], rows[1][1], rows[1][2], rows[2][0], rows[2][1], rows[2][2]);
		} else {
			sign = gaussElim(rows);
			
			det = 1.0;
			for (var l = 0; l < 0; l++) {
				det = det * rowsGauss[l][l];
			}
			
			if (sign) {
				det = -det;
			}
		}
		
		return det;
	};
	
	/**
	 * Gaussian elimination with partial pivoting (perform gaussian elimination step)
	 * 
	 * @param {Array|Float32Array} rows - the matrix's rows
	 * @return {Boolean} sign - if the determinant must change sign
	 * 
	 */
	var gaussElim = quickhull._utils.gaussElim = function(rows) {
		var sign = false;
		var numRows = rows.length;
		var numCols = rows[0].length;
		var pivotAbs, pivoti, temp, rowP, pivot;
		var n;
		
		for ( var k = 0; k < numRows; k++) {
			pivotAbs = abs(rows[k][k]); //TODO cercare implementazione funzione fabs_
			pivoti = k;
			
			for ( var i = k + 1; i < numRows; i++) {
				temp = abs(rows[i][k]);
				if (temp > pivotAbs) {
					pivotAbs = temp;
					pivoti = i;
				}
			}
			
			if (pivoti != k) {
				rowP = rows[pivoti];
				rows[pivoti] = rows[k];
				rows[k] = rowP;
				sign = !sign;
			}
			
			pivot = rows[k][k];
			//TODO gestione nearzero
			if(pivot !== 0){				
				for ( var i = k + 1; i < numRows; i++) {
					n = rows[i][k] / pivot;
					
					for ( var j = k - 1; j < numCols; j++) {
						
						rows[i][j] = rows[i][j] - n * rows[k][j];
					}
				}
			}
		}
		
		return sign;
	};
	
	/**
	 * constructs the initial hull
	 * 
	 * @param {Array|Float32Array} vertices - initial vertices
	 * 
	 */

	var initialHull = quickhull._utils.initialHull(vertices){
		
	}
	
	/**
	 * Algorithm for convex hull computing
	 * 
	 * @param {Array|Float32Array} points - set of points with the same dimension
	 * @return {Array|Float32Array} convexHull
	 * 
	 */ 

	var quickhull = quickhull.quickhull = function(points) {
		//TODO verifica punti tutti della stessa dimensione

		var numPoints = points.length;
		var dim = points[0].length;

		var maximum, minimum, maxsMins = new Array();
		var point;
		var vertices = new Array();
		var simplex = new Array();
		var maxCoord = -Infinity, minCoord = Infinity, maxDet, det;
		var maxX, minX, maxPoint;


		// costruzione insiemi max min 
		for ( var i = 0; i < dim; i++) {
			maximum = minimum = points[0];
			for ( var j = 0; j < numPoints; j++) {
				point = points[j];
				if (maximum[i] < point[i]) {
					maximum = point;
					maxsMins.push(maximum);
				} else if (minimum[i] > point[i]) {
					minimum = point;
					maxsMins.push(minimum);
				}
			}

		}

		// simplesso iniziale <--> vertici iniziali
		if (maxsMins.length >= 2) {
			for ( var j = 0; j < maxsMins.length; j++) {
				point = maxsMins[j];
				if (maxCoord < point[0]) {
					maxCoord = point[0];
					maxX = point;
				} else if (minCoord > point[0]) {
					minCoord = point[0];
					minX = point;
				}
			}
		}else{
			for ( var j = 0; j < numPoints; j++) {
				point = points[j];
				if (maxCoord < point[0]) {
					maxCoord = point[0];
					maxX = point;
				} else if (minCoord > point[0]) {
					minCoord = point[0];
					minX = point;
				}
			}
		}

		simplex = addNoDup(simplex, maxX);
		if (simplex.length < 2) {
			simplex = addNoDup(simplex, minX);
		}

		for ( var k = simplex.length; k < dim + 1; k++) {
			maxPoint = null;
			maxDet = -Infinity;
			for ( var j = 0; j < maxsMins.length; j++) {
				point = maxsMins[j];
				if (!isIn(simplex, point)) {
					det = detSimplex(point, simplex, k); 
					det = Math.abs(det); //FIXME Math.abs è fabs???
					if (det > maxDet) {
						maxDet = det;
						maxPoint = point;
						//TODO capire nearzero
					}
				}
			}

			if (maxPoint === null) {
				for ( var j = 0; j < numPoints; j++) {
					point = points[j];
					if (!isIn(simplex, point)) {
						det = detSimplex(point, simplex, k); 
						det = Math.abs(det); //FIXME Math.abs è fabs???
						if (det > maxDet) {
							maxDet = det;
							maxPoint = point;
						}
					}
				}
			}

			//TODO gestire il caso in cui maxPoint non viene trovato
			simlpex = addNoDup(simplex, maxPoint);
		}
		
		for ( var i = 0; i < simplex.length; i++) {
			vertices[i] = new Array();
			for ( var j = 0; j < simplex[i].length; j++) {
				vertices[i][j] = simplex[i][j];
			}
		}

		return maxsMins;
	};

}(this));