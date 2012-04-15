/*!
 * quickhull.js
 * Quickhull algoritm in javascript translated from 
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
	var sqrt = Math.sqrt;
	
	var facetId = 0.0;
	
	function Facet(){
		this.id = facetId++;
		this.neighbors = new Array();
		this.furthestDist = 0.0;
		this.bestPoint = new Array();
		this.simplicial = true;
		this.good = true;
		this.newfacet = true;
		this.vertices = new Array();
		this.topoOriented = true;
		this.normal = new Array();
		this.offset = 0.0;
		this.outsideSet = new Array();
	};

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
			
			equal = comparePoints(tempPoint, point);

			if (equal) {
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
	};
	
	/**
	 * Determinant of 2x2 matrix
	 * 
	 * @param {Number} a1 - coordinate 0,0
	 * @param {Number} a2 - coordinate 0,1
	 * @param {Number} b1 - coordinate 1,0
	 * @param {Number} b2 - coordinate 1,1
	 * @return {Number} determinant
	 * 
	 */
	
	var det2 = quickhull._utils.det2 = function(a1, a2, b1, b2) {
		return (a1 * b2) - (a2 * b1);
	};
	
	/**
	 * Determinant of 3x3 matrix
	 * 
	 * @param {Number} a1 - coordinate 0,0
	 * @param {Number} a2 - coordinate 0,1
	 * @param {Number} a2 - coordinate 0,2
	 * @param {Number} b1 - coordinate 1,0
	 * @param {Number} b2 - coordinate 1,1
	 * @param {Number} b2 - coordinate 1,2
	 * @param {Number} c1 - coordinate 2,0
	 * @param {Number} c2 - coordinate 2,1
	 * @param {Number} c2 - coordinate 2,2
	 * @return {Number} determinant
	 * 
	 */
	
	var det3 = quickhull._utils.det3 = function(a1, a2, a3, b1, b2, b3, c1, c2, c3) {
		return (a1) * det2(b2,b3,c2,c3) - (b1) * det2(a2,a3,c2,c3) + (c1)*det2(a2,a3,b2,b3);
	};

	/**
	 * Compute the determinant of a simplex
	 * 
	 * @param {Array|Float32} point - point apex
	 * @param {Array|Float32Array} simplex - base points
	 * @param {Number} dim - dimension
	 */

	var detSimplex = quickhull._utils.detSimplex = function(point, simplex, dim) {
		
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
			sign = gaussElim(rows, false);
			
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
	var gaussElim = quickhull._utils.gaussElim = function(rows, sign) {
		//var sign = false;
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
	 * @return {Array|Facet} initial hull
	 * 
	 */

	var initialHull = quickhull._utils.initialHull = function (vertices) {
		var topoOrient = true;
		var newFacets = new Array();
		var newFacet;
		var interiorPoint;
		var dist;
		
		for ( var i = 0; i < vertices.length; i++) {
			newFacet = new Facet();
			newFacet.topoOriented = topoOrient;
			//TODO ordinare i vertici
			for ( var j = 0; j < vertices.length; j++) {
				if (j != i) {
					newFacet.vertices.push(vertices[j]);
				}
			}
			newFacets.push(newFacet);
			topoOrient = !topoOrient;
		}
		
		for ( var i = 0; i < newFacets.length; i++) {
			newFacet = newFacets[i];
			for ( var j = 0; j < newFacets.length; j++) {
				if (newFacet.id !== newFacets[j].id) {
					newFacet.neighbors.push(newFacets[j]);
				}
			}
		}
		
		/*console.log("Facce: ");
		for ( var i = 0; i < newFacets.length; i++) {
			console.log("faccia.id = " + newFacets[i].id);
			console.log("faccia.topoOriented = " + newFacets[i].topoOriented);
			console.log("faccia.vertices: ");
			for ( var j = 0; j < newFacets[i].vertices.length; j++) {
				console.log("[" + newFacets[i].vertices[j] + "]");	
			}
			var testo = "";
			for ( var j = 0; j < newFacets[i].neighbors.length; j++) {
				if (j === 0) {
					testo += "" + newFacets[i].neighbors[j].id;
				} else {
					testo += "," + newFacets[i].neighbors[j].id;
				}
			}
			console.log("faccia.neighbors: [" + testo + "]");
		}*/
		
		//TODO verificare reset list
		interiorPoint = getCenter(vertices);
		console.log("Center: [" + interiorPoint + "]");
		
		setFacetPlane(newFacets[0]);
		
		dist = distPlane(interiorPoint, newFacets[0]);
		
		if (dist > 0) {
			for ( var i = 0; i < newFacets.length; i++) {
				newFacets[i].topoOriented = !(newFacets[i].topoOriented);
			}
		}
		
		for ( var i = 0; i < newFacets.length; i++) {
			setFacetPlane(newFacets[i]);
		}
		//TODO implementare checkflipped, gestione minangle, qh_checkpolygon, qh_checkconvex
		
		return newFacets;
	};
	
	/**
	 * return distance from point to facet
	 * 
	 * @param {Array|Float32Array} point - point
	 * @param {Facet} facet - facet
	 * @return {Number} distance of point from facet
	 */

	var distPlane = quickhull._utils.distPlane = function (point, facet) {
		var dist = facet.offset, normal = facet.normal;
		
		for ( var i = 0; i < normal.length; i++) {
			dist = dist + (point[i]*normal[i]);
		}
		
		return dist;
	};
	
	/**
	 * returns arithmetic center of a set of vertices
	 * 
	 * @param {Array|Float32Array} vertices - set of vertices
	 * @return {Float32Array} center point
	 */

	var getCenter = quickhull._utils.getCenter = function (vertices) {
		
		var center = new Array();
		var size = vertices.length;
		var dim = vertices[0].length;
		var coordinata;
		
		//TODO implementare controllo che siano almeno due punti
		for ( var i = 0; i < dim; i++) {
			coordinata = 0.0;
			for ( var j = 0; j < vertices.length; j++) {
				coordinata = coordinata + vertices[j][i];
			}
			center.push(coordinata/size);
		}
		
		return center;
	};
	
	/**
	 * sets the hyperplane for a facet and the offset
	 * 
	 * @param {Facet} facet - the facet to set hyperplane
	 * 
	 */

	var setFacetPlane = quickhull._utils.setFacetPlane = function (facet) {
		var point0 = facet.vertices[0], normal = new Array();
		var offset, rows = new Array();
		
		if (point0.lenght <= 4) {
			for ( var i = 0; i < facet.vertices.length; i++) {
				rows.push(facet.vertices[i]);
			}
			
			normal = hyperplaneDet(rows, point0, facet.topoOriented);
		}
		
		if (point0.lenght > 4) { //TODO || nearzero
			rows = new Array();
			//TODO gestire ciclo:
			//FOREACHvertex_(facet->vertices) {
		    //  if (vertex->point != point0) {
		    //      qh gm_row[i++]= gmcoord;
			//      coord= vertex->point;
			//      point= point0;
			//      for (k=qh hull_dim; k--; )
			//        *(gmcoord++)= *coord++ - *point++;
			//    }
			//}
			for ( var i = 0; i < facet.vertices.length; i++) {
				rows.push(facet.vertices[i]);
			}
			
			normal = hyperplaneGauss(rows, point0, facet.topoOriented);
		}
		
		facet.normal = normal;
		offset = -(point0[0] * normal[0]);
		
		for ( var i = 1; i < normal.length; i++) {
			offset = offset - (point0[i] * normal[i]);
		}
		
		facet.offset = offset;
	};
	
	/**
	 * given an upper-triangular rows array and a sign,
     * solve for normal equation x using back substitution over rows U
	 * 
	 * @param {Array|Float32Array} rows - one row per point
	 * @param {Number} numrow - number of rows
	 * @param {Number} numcol - number of cols
	 * @param {Boolean} sign - 
	 * @return {Float32Array} normal normalized
	 * 
	 */

	var backNormal = quickhull._utils.backNormal = function (rows, numrow, numcol, sign) {
		var normal = new Array(numcol), diagonal;
		normal[normal.length - 1] = (sign ? -1.0 : 1.0);
		for ( var i = numrow - 2; i >= 0; i--) {
			normal[i] = 0.0;
			for ( var j = i + 1; j < numcol; j++) {
				normal[i] = normal[i] - (rows[i][j] * normal[j]);
			}
			diagonal = rows[i][i];
			//TODO gestire caso diagonal uguale a zero
			normal[i] = normal[i]/diagonal;
		}
		
		return normal;
	};
	
	/**
	 * return normalized hyperplane equation from oriented simplex
	 * 
	 * @param {Array|Float32Array} rows - one row per point
	 * @param {Float32Array} point0 - any row
	 * @param {Boolean} topoOrient - flips all signs
	 * @return {Float32Array} normal normalized
	 * 
	 */

	var hyperplaneGauss = quickhull._utils.hyperplaneGauss = function (rows, point0, topoOrient) {
		var dim = rows[0].length;
		var normal;
		var sign = gaussElim(rows, topoOrient);
		
		for ( var i = 0; i < dim; i++) {
			if (rows[i][i] < 0) {
				sign = !sign;
			}
		}
		
		normal = backNormal(rows, rows.length - 1, rows[0].length, sign);
		//TODO gestione nearzero
		
		normalize2(normal, dim, topoOrient);
		
		return normal;
	};
	
	/**
	 * return normalized hyperplane equation from oriented simplex
	 * 
	 * @param {Array|Float32Array} rows - one row per point
	 * @param {Float32Array} point0 - any row
	 * @param {Boolean} topoOrient - flips all signs
	 * @return {Float32Array} normal normalized
	 * 
	 */

	var hyperplaneDet = quickhull._utils.hyperplaneDet = function (rows, point0, topoOrient) {
		/*given two indices into rows[],

         compute the difference between X, Y, or Z coordinates
		  #define dX( p1,p2 )  ( *( rows[p1] ) - *( rows[p2] ))
		  #define dY( p1,p2 )  ( *( rows[p1]+1 ) - *( rows[p2]+1 ))
		  #define dZ( p1,p2 )  ( *( rows[p1]+2 ) - *( rows[p2]+2 ))
		  #define dW( p1,p2 )  ( *( rows[p1]+3 ) - *( rows[p2]+3 ))*/
		var dim = rows[0].length;
		var normal = new Array();
		
		if (dim === 2) {
			normal.push(rows[1][1] - rows[0][1]);//dY(1,0)
			normal.push(rows[0][0] - rows[1][0]);//dX(0,1)
			normalize2(normal, dim, topoOrient);
		} else if (dim === 3) {
			normal.push(det2(rows[2][1] - rows[0][1],//dY(2,0) 
					         rows[2][2] - rows[0][2],//dZ(2,0)
					         rows[1][1] - rows[0][1],//dY(1,0)
					         rows[1][2] - rows[0][2]));//dZ(1,0)
			
			normal.push(det2(rows[1][0] - rows[0][0],//dX(1,0)
			         		 rows[1][2] - rows[0][2],//dZ(1,0)
			         		 rows[2][0] - rows[0][0],//dX(2,0)
			         		 rows[2][2] - rows[0][2]));//dZ(2,0)
			
			normal.push(det2(rows[2][0] - rows[0][0],//dX(2,0)
			         		 rows[2][1] - rows[0][1],//dY(2,0)
			         		 rows[1][0] - rows[0][1],//dX(1,0)
			         		 rows[1][1] - rows[0][1]));//dY(1,0)
			
			normalize2(normal, dim, topoOrient);
		} else if (dim === 4) {
			normal.push(-det3(rows[2][1] - rows[0][1],//dY(2,0)
	         		 		  rows[2][2] - rows[0][2],//dZ(2,0)
	         		 		  rows[2][3] - rows[0][3],//dW(2,0)
	         		 		  rows[1][1] - rows[0][1],//dY(1,0)
	         		 		  rows[1][2] - rows[0][2],//dZ(1,0)
	         		 		  rows[1][3] - rows[0][3],//dW(1,0)
	         		 		  rows[3][1] - rows[0][1],//dY(3,0)
	         		 		  rows[3][2] - rows[0][2],//dZ(3,0)
	         		 		  rows[3][3] - rows[0][3]));//dW(3,0)
			
			normal.push(det3(rows[2][0] - rows[0][0],//dX(2,0) 
   		 		  			  rows[2][2] - rows[0][2],//dZ(2,0)
   		 		  			  rows[2][3] - rows[0][3],//dW(2,0)
   		 		  			  rows[1][0] - rows[0][0],//dX(1,0) 
   		 		  			  rows[1][2] - rows[0][2],//dZ(1,0)
   		 		  			  rows[1][3] - rows[0][3],//dW(1,0)
   		 		  			  rows[3][0] - rows[0][0],//dX(3,0) 
   		 		  			  rows[3][2] - rows[0][2],//dZ(3,0)
   		 		  			  rows[3][3] - rows[0][3]));//dW(3,0)
			
			normal.push(-det3(rows[2][0] - rows[0][0],//dX(2,0) 
   		 		  			  rows[2][1] - rows[0][1],//dY(2,0)
   		 		  			  rows[2][3] - rows[0][3],//dW(2,0)
   		 		  			  rows[1][0] - rows[0][0],//dX(1,0) 
   		 		  			  rows[1][1] - rows[0][1],//dY(1,0)
   		 		  			  rows[1][3] - rows[0][3],//dW(1,0)
   		 		  			  rows[3][0] - rows[0][0],//dX(3,0) 
   		 		  			  rows[3][1] - rows[0][1],//dY(3,0)
   		 		  			  rows[3][3] - rows[0][3]));//dW(3,0)
			
			normal.push(det3(rows[2][0] - rows[0][0],//dX(2,0)
	 		  			  	 rows[2][1] - rows[0][1],//dY(2,0)
	 		  			  	 rows[2][2] - rows[0][2],//dZ(2,0)
	 		  			  	 rows[1][0] - rows[0][0],//dX(1,0)
	 		  			  	 rows[1][1] - rows[0][1],//dY(1,0)
	 		  			  	 rows[1][2] - rows[0][2],//dZ(1,0)
	 		  			  	 rows[3][0] - rows[0][0],//dX(3,0)
	 		  			  	 rows[3][1] - rows[0][1],//dY(3,0)
	 		  			  	 rows[3][2] - rows[0][2]));//dZ(3,0)
			
			normalize2(normal, dim, topoOrient);
		}
		return normal;
	};
	
	/**
	 * normalize a vector
	 * 
	 * @param {Float32Array} normal - vector to normalize
	 * @param {Boolean} topoOrient - flips all signs
	 * 
	 */
	
	var normalize2 = quickhull._utils.normalize2 = function (normal, dim, topoOrient) {
		var somma = 0.0, norm = 0.0, temp;
		for ( var i = 0; i < dim; i++) {
			somma = somma + (normal[i] * normal[i]);
		}
		
		norm = sqrt(somma);
		
		if (norm === 0.0) {
			temp = sqrt(1.0/dim);
			for ( var i = 0; i < dim; i++) {
				normal[i] = temp;
			}
		} else {
			if (!topoOrient) {
				norm = -norm;
			}
			for ( var i = 0; i < dim; i++) {
				//TODO gestire divisioni per zero
				normal[i] = normal[i] / norm;
			}
		}
	};
	
	/**
	 * Compare points with same dimension
	 * 
	 * @param {Float32Array} point1 - first point
	 * @param {Float32Array} point2 - second point
	 * @return {Boolean} if equal
	 * 
	 */ 

	var comparePoints = quickhull._utils.comparePoints = function(point1, point2) {
		var equal = true;
		
		for ( var j = 0; j < point1.length && equal; j++) {
			// se almeno una coordiata è differente non sono lo stesso punto
			if (point1[j] !== point2[j]) {
				equal = false;
			}
		}
		
		return equal;
	};
	
	/**
	 * partitions all points in points/numpoints to the outsidesets of facets
	 * 
	 * @param {Array|Facet} facets - set of facets
	 * @param {Array|Float32Array} vertices - set of points facets' vertices
	 * @param {Array|Float32Array} points - set of points
	 * @param {Number} - number of points
	 * 
	 */ 
	
	var partitionAll = quickhull._utils.partitionAll = function(facets, vertices, points, numPoints) {
		var pointSet = new Array();
		var isVertex;
		var bestPoint, bestDist, dist;
		
		for ( var i = 0; i < points.length; i++) {
			isVertex = isIn(vertices, points[i]);
			if (!isVertex) {
				pointSet.push(points[i]);
			}
		}
		
		for ( var i = 0; i < facets.length; i++) {
			bestPoint = null; 
			bestDist = 0.0;
			
			for ( var j = 0; j < pointSet.length; j++) {
				dist = distPlane(pointSet[j], facets[i]);
				if (dist > 0) {//TODO verificare le condizioni di punto esterno alla faccia
					
					if (bestPoint === null) {
						bestPoint = pointSet[j];
						bestDist = dist;
					} else if (dist > bestDist) {
						facets[i].outsideSet.push(bestPoint);
						bestPoint = pointSet[j];
						bestDist = dist;
					} else {
						facets[i].outsideSet.push(pointSet[j]);
					}
				} /*else if (dist === 0) {
				//TODO testare profondamente metodo .splice
					pointSet.splice(i, 1);
				}*/
				if (bestPoint !== null) {
					facets[i].furthestDist = bestDist;
					facets[i].bestPoint = bestPoint;
				}
			}
		}
	};
	
	/**
	 * return the facet with furthest of all furthest points searches all facets
	 * 
	 * @param {Array|Facet} facets - set of facets
	 * @return {Facet} the facet with furthest 
	 * 
	 */ 
	
	var furthestNext = quickhull._utils.furthestNext = function(facets) {
		var bestDist = bestFacet.furthestDist, bestFacet = facets[0];
		
		for ( var i = 1; i < facets.length; i++) {
			if (bestDist < facets[i].furthestDist) {
				bestDist = facets[i].furthestDist;
				bestFacet = facets[i];
			}
		}
		
		return bestFacet;
		/*facetT *facet, *bestfacet= NULL;
		  realT dist, bestdist= -REALmax;

		  FORALLfacets {
		    if (facet->outsideset) {
		#if qh_COMPUTEfurthest
		      pointT *furthest;
		      furthest= (pointT*)qh_setlast(facet->outsideset);
		      zinc_(Zcomputefurthest);
		      qh_distplane(furthest, facet, &dist);
		#else
		      dist= facet->furthestdist;
		#endif
		      if (dist > bestdist) {
		        bestfacet= facet;
		        bestdist= dist;
		      }
		    }
		  }
		  if (bestfacet) {
		    qh_removefacet(bestfacet);
		    qh_prependfacet(bestfacet, &qh facet_next);
		    trace1((qh ferr, 1029, "qh_furthestnext: made f%d next facet(dist %.2g)\n",
		            bestfacet->id, bestdist));
		  }*/
	};
	
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
		var convexHull = new Array();
		var maxCoord = -Infinity, minCoord = Infinity, maxDet, det;
		var maxX, minX, maxPoint;
		var initHull;
		var nextFacet;

		console.log("Inizio calcolo convex hull con i punti: " + points);

		// costruzione insiemi max min per ogni dimesione
		for ( var i = 0; i < dim; i++) {
			maximum = minimum = points[0];
			for ( var j = 0; j < numPoints; j++) {
				point = points[j];
				if (maximum[i] < point[i]) {
					maximum = point;
				} else if (minimum[i] > point[i]) {
					minimum = point;
				}
			}
			maxsMins.push(maximum);
			maxsMins.push(minimum);
		}
		
		/*console.log("Insieme iniziale Maxs-Mins: ");
		for ( var i = 0; i < maxsMins.length; i++) {
			console.log("[" + maxsMins[i] + "]");
		}*/

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

		addNoDup(simplex, maxX);
		//console.log("maxX: " + maxX);
		//console.log("minX: " + minX);
		if (simplex.length < 2) {
			addNoDup(simplex, minX);
		}
		
		/*console.log("Simplesso primi punti: ");
		for ( var i = 0; i < simplex.length; i++) {
			console.log("[" + simplex[i] + "]");
		}*/

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
		
		/*console.log("Simplesso/vertici iniziale/i: ");
		for ( var i = 0; i < simplex.length; i++) {
			console.log("[" + simplex[i] + "]");
		}*/
		
		for ( var i = 0; i < simplex.length; i++) {
			vertices[i] = new Array();
			for ( var j = 0; j < simplex[i].length; j++) {
				vertices[i][j] = simplex[i][j];
			}
		}
		
		initHull = initialHull(vertices);
		
		partitionAll(vertices, points, numPoints);
		
		nextFacet = furthestNext(initHull);

		return convexHull;
	};

}(this));