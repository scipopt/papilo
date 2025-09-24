import re

def parse_logfile(filename1, filename2):
	results = []
	problem_re = re.compile(r"read problem </scratch/htc/jholly/scip/check/IP/miplib2017-benchmark/(.+?)\.mps\.gz>")
	presolver_re = re.compile(r"\(\s*([0-9]*\.?[0-9]+)s\)\s*MILP presolver", re.IGNORECASE)
	nothing_re = re.compile(r"found nothing")
	changes_re = re.compile(r"(\d+) aggregations, (\d+) fixings, (\d+) bound changes")

	with open(filename1, "r") as f:
		lines1 = f.readlines()
	with open(filename2, "r") as g:
		lines2 = g.readlines()	

	i = 0
	while i < len(lines1):
		line = lines1[i]
		match_problem = problem_re.search(line)
		if match_problem:
			X = match_problem.group(1)
			i += 1
			while i < len(lines1) and not presolver_re.search(lines1[i]):
				i += 1
			if i < len(lines1):
				match_presolver = presolver_re.search(lines1[i])
				if match_presolver:
					Y = float(match_presolver.group(1))

					if i < len(lines1) and nothing_re.search(lines1[i]):
						results.append([X, Y, 0, 0, 0])
					elif i < len(lines1):
						match_changes = changes_re.search(lines1[i])
						if match_changes:
							z, v, w = map(int, match_changes.groups())
							results.append([X, Y, v, z, w])
		i += 1
	i = 0
	instance = 0
	while i < len(lines2) and instance < len(results):
		line = lines2[i]
		match_problem = problem_re.search(line)
		if match_problem and match_problem.group(1) == results[instance][0] :
			X = match_problem.group(1)
			i += 1
			while i < len(lines2) and not presolver_re.search(lines2[i]):
				i += 1

			if i < len(lines2):
				match_presolver = presolver_re.search(lines2[i])
				if match_presolver:
					Y = float(match_presolver.group(1))

					if i < len(lines2) and nothing_re.search(lines2[i]) :
						results[instance] += [Y, 0, 0, 0]
						i=0
						instance+=1
					elif i < len(lines2) :
						match_changes = changes_re.search(lines2[i])
						if match_changes:
							z, v, w = map(int, match_changes.groups())
							results[instance] += [Y, v, z, w]
							instance+=1
							i=0
							
		i += 1

	i = 0
	outputstring = ""
	while i < len(results) :
		if results[i][2] == results[i][6] and results[i][3] == results[i][7] and results[i][4] == results[i][8] and ( results[i][1] + 0.0001 ) /( results[i][5] + 0.0001 ) < 1.5 and ( results[i][1] + 0.0001 ) / (results[i][5] + 0.0001) > 0.66 :
			i += 1
			continue
		outputstring += results[i][0]
		outputstring += " & "
		outputstring += str(results[i][1])
		outputstring += " & "
		outputstring += str(results[i][2])
		outputstring += " & "
		outputstring += str(results[i][3])
		outputstring += " & "
		outputstring += str(results[i][4])
		outputstring += " & "
		outputstring += str(results[i][5])
		outputstring += " & "
		outputstring += str(results[i][6])
		outputstring += " & "
		outputstring += str(results[i][7])
		outputstring += " & "
		outputstring += str(results[i][8])
		outputstring += " \\\\ \n"
		i += 1
	print( outputstring )
	return

parse_logfile(r"check.miplib2017_benchmark.scip.M640.default-s3.out", r"check.miplib2017_benchmark.scip.M640.clique-s3.out")