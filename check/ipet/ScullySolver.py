import re

from ipet.parsing.Solver import SCIPSolver

VALUE = "value"
RUNTIME = "runtime"
FIRST_SOLUTION = "first"
BEST_SOLUTION = "best"

DEFAULT_VALUE = -1.0

class ScullySolver(SCIPSolver):
    solverId = "Scully"
    # used to identify the solver
    recognition_expr = re.compile("Scully version")
    version_expr = re.compile("Scully version (\S+)")

    floating_point_expr = "[-+]?[0-9]*\.?[0-9]*"

    value_expr = re.compile("Algorithm with objective value\s+(\S+)")
    runtime_expr = re.compile("Algorithm with objective value\s+" + floating_point_expr + "\s+finished after\s+(\S+)")
    first_solution_expr = re.compile("First solution found in round\s+(\S+)")
    best_solution_expr = re.compile("Best solution found in round\s+(\S+)")

    def __init__(self, **kw):
        super(ScullySolver, self).__init__(**kw)

    def reset(self):
        super(ScullySolver, self).reset()
        self.addData(VALUE, DEFAULT_VALUE)
        self.addData(RUNTIME, DEFAULT_VALUE)
        self.addData(FIRST_SOLUTION, DEFAULT_VALUE)
        self.addData(BEST_SOLUTION, DEFAULT_VALUE)

    def extractPrimalboundHistory(self, line):
        pass

    def extractDualboundHistory(self, line):
        pass

    def extractHistory(self, line):
        pass

    def extractOptionalInformation(self, line: str):
        self.extractByExpression(line, self.value_expr, VALUE)
        self.extractByExpression(line, self.runtime_expr, RUNTIME)
        self.extractByExpression(line, self.best_solution_expr, BEST_SOLUTION)
        self.extractByExpression(line, self.first_solution_expr, FIRST_SOLUTION)
