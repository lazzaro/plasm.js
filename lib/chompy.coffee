{PI, E, log, sin, cos, tan, asin, acos, atan, atan2, ceil, floor, sqrt, exp, abs, round} = Math
apply = (args) ->  [f, x] = args; f.apply(null,[x])
comp2 = (f,g) -> (x) -> f (g x)
comp = (flist) -> flist.reduceRight comp2
cons = (flist) -> (x) -> flist.map (f) -> f x
cat = (a) -> [].concat a...
id = (a) -> a
k = (a) -> (b) -> a
aa = (f) -> (list) -> list.map (e) -> f e
distr = (args) -> [list,e] = args; [el,e] for el in list
distl = (args) -> [e,list] = args; [e,el] for el in list
insr = (f) -> (seq) -> seq.reduceRight f
insl = (f) -> (seq) -> seq.reduce f
bigger = (a,b) -> if a > b then a else b
smaller = (a,b) -> if a < b then a else b
biggest = (args) -> args.reduce bigger #see below
biggest = (args) -> (insr bigger) args #see above
smallest = (args) -> (insr smaller) args 
list = (args) -> (cons [id]) args
len = (args) -> args.length
reverse = (args) -> if args.length > 1 then (args[i] for i in [args.length-1..0]) else args
tail = (args) -> if args.length > 0 then args.splice 1, args.length-1 else args
butlast = (args) -> if args.length > 1 then reverse tail reverse args else []
al = (args) -> cat args
ar = (args) -> cat args
repeat = (n) -> (args) -> (args for i in [0...n])
replica = (n) -> (args) -> cat (args for i in [0...n])
sum = (args) -> if typeof args[0] is 'number' then (insl (x,y) -> x+y) args else aa(insl (x,y) -> x+y)(trans args)
sub = (args) -> if typeof args[0] is 'number' then (insl (x,y) -> x-y) args else aa(insl (x,y) -> x-y)(trans args)
mul = (args) -> if typeof args[0] is 'number' then (insl (x,y) -> x*y) args else aa(insl (x,y) -> x*y)(trans args)
div = (args) -> if typeof args[0] is 'number' then (insl (x,y) -> x/y) args else aa(insl (x,y) -> x/y)(trans args)
trans = (args) -> 
	n = args.length; m = args[0].length; args = cat args
	((args[j*m+i] for j in [0...n]) for i in [0...m])
vect = (binaryop) -> (args) -> aa(binaryop) trans args
myprint = (string,params) -> console.log string, params, "\n"
mat = (m,n) -> (args) -> ((args[i*n+j] for j in [0...n]) for i in [0...m])
isNumber = (n) -> (not isNaN parseFloat n) and isFinite n

#///////////////////////////////////////////////////////////////////////////////

progressive_sum = (args) -> 
	al [0, (insr (x,y) -> x+y) args[0..i] for i in [0...args.length]]

#///////////////////////////////////////////////////////////////////////////////

type = (obj) ->
	if obj == undefined or obj == null
		return String obj
	classToType = new Object
	for name in "Boolean Number String Function Array Date RegExp".split(" ")
		classToType["[object " + name + "]"] = name.toLowerCase()
	myClass = Object.prototype.toString.call obj
	if myClass of classToType
		return classToType[myClass]
	return "object"
typedPrint = (args) -> 
	console.log "#{type args}::#{args}";
	args
clone = (obj) ->
	if not obj? or typeof obj isnt 'object' 
		return obj
	newInstance = new obj.constructor()
	for key of obj
		newInstance[key] = clone obj[key]
	newInstance

#///////////////////////////////////////////////////////////////////////////////

PRECISION = 1E7
fixedPrecision = (number) -> 
	int = (if number>0 then floor else ceil) number
	number = (if number>0 then ceil else floor)(PRECISION * number) / PRECISION
	if abs(number-int) <= 1.0/PRECISION then int else number

fcode = (point) -> (aa fixedPrecision) point
code = (point) -> "[#{fcode point}]"

decode = (string) -> +string  # => a number
uncode = (pointCode) -> (aa decode) pointCode.split(',')

string2numberList = (string) ->
	if string is '[]' then [] else
		regex = /[,|\[|\]]/ # regex => "[" or "," or "]"
		(aa Number) butlast tail string.split(regex)

mapcomp = (map1,map2) -> map = {}; map[k] = map2[v] for k,v of map1

#///////////////////////////////////////////////////////////////////////////////

revert = (cell) -> 
	len = cell.length
	if len >1 then cat [cell[len-1], cell[1...len-1], cell[0]] else cell
remove_duplicates  = (hasDupes) -> 
	dict = {}; (dict[code(item)] = item for item in hasDupes \
		when not dict[code(revert item)]? and not dict[code(item)]?)
rotate = (cell) -> 
	if cell.length > 1 then cat [cell[1...cell.length],[cell[0]]] else cell
facets  = (cell) -> 
	out = []; for h in [0...cell.length]
		facet = (k for k,i in cell when i isnt h)
		out.push  if h%2 is 1 then revert facet else facet
	out
skeleton  = (h_cells) -> 
	remove_duplicates cat (facets cell for cell in h_cells)
cell_complex = (d_cells) -> 
	if d_cells.length > 0  
		dim = d_cells[0].length-1 
		cells = new Array(dim)
		cells[dim] = d_cells
		cells[h-1] = skeleton cells[h] for h in [dim..1]
		cells
	else
		dim = -1
		cells = []
mkCellDB  = (complex) -> 
	complex = complex or []
	dictos = []
	for skel,d in complex
		dictos[d] = {}; dictos[d][code(cell)] = k for cell,k in skel
	dictos
homology_maps = (dictos) -> 
	if dictos.length > 0
		dim = dictos.length-1; d = 1
		homology = ([] for i in [0..dim])
		if dim > 1
			skel = dictos[1]
			for cell of skel
				simplex = string2numberList cell
				homology[1].push ([skel[cell], facet[0]] for facet in facets simplex)	
			homology[1] = cat homology[1]
			for skel in dictos[2..dim]
				d += 1; 
				for cell of skel
					for facet in facets string2numberList cell
						if dictos[d-1][code(facet)]?
							key = dictos[d-1][code(facet)] 
						else 
							key = dictos[d-1][code(revert(facet))]						
						homology[d].push [skel[cell], key]
		homology
	else []

#///////////////////////////////////////////////////////////////////////////////

coords_distribute = (x) -> 
	out = cat( aa(ar)(distr(e)) for e in x)

subcomplex = (d,args) ->
    (args[i...i+d] for i in [0..args.length-d])

shift = (n, listoflists) -> (x+n for x in seq) for seq in listoflists

#///////////////////////////////////////////////////////////////////////////////

class PointSet
	constructor: (points) ->
		points = points or []
		if points.length > 0
			@rn = points[0].length
			@dict = {}; @dict[fcode(point)] = i for point,i in points
			map1 = {}; map1[i] = decode(@dict[fcode(point)])  for point,i in points
			[map2,k] = [{},0]; for pcode,pid of @dict
				map2[decode(pid)] = k; k += 1
			@map = (v for k,v of mapcomp map1,map2)
			k = 0; for pcode of @dict
				@dict[pcode] = k; k += 1
			@verts = (uncode pcode for pcode of @dict)
			@m = @verts.length
		else
			@rn = 0
			@dict = {}
			@map = []
			@verts = []
			@map = 0
		
	update: (modify) ->
		@verts[pid] = modify (uncode pcode) for pcode,pid of @dict
		@dict = {}; (@dict[code(point)] = pid for point,pid in @verts)
		
	t: (indices,values) -> 
		vect = (0 for k in [0..@rn])
		vect[indices[h]] = values[h] for h in [0...indices.length]
		@update (point) -> sum([point, vect])
		@
		
	s: (indices,values) -> 
		vect = (1 for k in [0..@rn])
		vect[indices[h]] = values[h] for h in [0...indices.length]
		@update (point) -> mul([point, vect]) 
		@
		
	###
	r: (axes, angle) -> 
		@update (point) -> id([axes, angle]) 
		@
	###

class Topology
	constructor: (vertices,d_cells) ->
		vertices = vertices or []
		d_cells = d_cells or []
		@dim = if d_cells.length > 0 then d_cells[0].length-1 else -1
		d_cells = ((vertices.map[k] for k in cell) for cell in d_cells)
		@dictos = mkCellDB cell_complex d_cells
		@homology = homology_maps @dictos
		@cells = (string2numberList cell for cell of dict for dict in @dictos)


class SimplicialComplex
	constructor: (points,d_cells) -> 
		points = points or []
		d_cells = d_cells or []
		@vertices = new PointSet(points)
		@faces = new Topology(@vertices,d_cells)

	t: (indices,values) -> @vertices.t(indices,values); @
	s: (indices,values) -> @vertices.s(indices,values); @

	r: (axes, angle) ->
		@vertices.r(axes, angle); 
		@

	extrude: (hlist) ->
		cells = clone @faces.cells
		dim = clone @faces.dim
		verts = clone @vertices.verts
		lastcoords = progressive_sum (aa)(abs)(hlist)
		if dim is 0
			cells = [[],[]]
			vertices = (aa list) lastcoords
			cells[1] = ([i,i+1] for i in [0...hlist.length+1])
		else
			simplexes = cells[dim]
			nverts = verts.length
			nsteps = lastcoords.length
			sliced_vertices = (replica nsteps) [verts]			
			vertices = coords_distribute(trans([repeat(nsteps)(verts), lastcoords]))
			extruded_simplices = []
			for cell in simplexes
				vertPtrs = cat([cell, cell.map (x) -> x+nverts])	 
				extruded_simplices.push subcomplex(dim+2,vertPtrs)
			final_simplices = []
			for i in [0..nsteps] 
				if hlist[i] > 0
					simplex_layer = shift nverts*i, cat extruded_simplices
					final_simplices.push simplex_layer
			cells = cat final_simplices
			#console.log  cells
		new SimplicialComplex(vertices, cells)

#///////////////////////////////////////////////////////////////////////////////


t = (indices,values) -> (obj) -> (clone obj).t(indices,values)
s = (indices,values) -> (obj) -> (clone obj).s(indices,values)


#///////////////////////////////////////////////////////////////////////////////


centroid = (obj) -> (face) ->
	A = (obj.vertices.verts[v]  for v in face)
	C = repeat(A.length)(1.0/A.length)
	point = numeric.dot(C,A)

###
obj = new SimplicialComplex [[0,0,0],[1,0,0],[0,1,0],[0,0,1]],[[0,1,2,3]]
###

#///////////////////////////////////////////////////////////////////////////////

simplexGrid = (args) ->
	hlist = args[0]
	lastcoords = progressive_sum aa(abs)(hlist)
	verts = aa(list)(lastcoords)
	cells = ([i,i+1] for i in [0...hlist.length] when hlist[i] > 0 )
	complex = new SimplicialComplex(verts,cells)
	for hlist in args[1...args.length]
		complex = complex.extrude(hlist)
	complex
	

#///////////////////////////////////////////////////////////////////////////////

free = (obj) -> 
	d = obj.faces.dim
	simplices = (obj.vertices.verts[k] for k in cell for cell in obj.faces.cells[d])
	out = []; for simplex in simplices
		outsimplex = new SimplicialComplex(simplex,[[0..d]])
		out.push outsimplex
	out

#///////////////////////////////////////////////////////////////////////////////

explode = (args) -> (scene) ->
	face = () -> item.faces.cells[item.faces.dim][0]
	centers = (centroid(item)(face()) for item in scene)
	scaledCenters = (mul([args,center]) for center in centers)
	translVectors = (sub(pair) for pair in trans([scaledCenters, centers]))
	scene[k].t([0...v.length],v) for v,k in translVectors

#///////////////////////////////////////////////////////////////////////////////

outline = (dim) -> (pol) ->
	#console.log pol
	verts = pol.vertices.verts
	faces_d = pol.faces.cells[dim]
	#console.log "verts",verts
	#console.log "faces_d",faces_d
	new SimplicialComplex(verts,faces_d)

#///////////////////////////////////////////////////////////////////////////////


#///////////////////////////////////////////////////////////////////////////////

obj = simplexGrid ([[1],[1],[1]]) 
obj = outline(2)(obj) ## OK


obj = outline(1)(obj) ## KO

obj = simplexGrid ([[1],[1]]) 
obj = outline(1)(obj) ## KO

console.log "outverts",obj.vertices.verts
console.log "outfaces",obj.faces.cells[obj.faces.dim]
model = viewer.draw(obj) 

