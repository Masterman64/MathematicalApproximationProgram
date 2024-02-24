# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 14:35:33 2023

@author: jones
"""

from sympy import sympify
from sympy import lambdify
from sympy import Symbol
from sympy import SympifyError
from sympy import solve
from sympy import nsolve

class StandardFunctionResult:
    def __init__(self, functionString, initialGuess, realRoot):
        self.functionString = functionString
        self.function = sympify(functionString)
        self.initialGuess = initialGuess
        self.realRoot = realRoot
        
        self.iterations = 0
        self.approximateRoot = 0.0
        self.forwardError = 0.0
        self.backwardsError = 0.0
        self.convergenceRate = 0.0
        
class BisectionMethodFunctionResult:
    def __init__(self, functionString, initialGuess, lowerBound, upperBound):
        self.functionString = functionString
        self.function = sympify(functionString)
        self.initialGuess = initialGuess
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        
        self.iterations = 0
        self.approximateRoot = 0.0
        self.intervalLength = 0.0
        self.error = 0.0
        self.convergenceRate = 0.0
    
def bisectionMethod(variable, function, lowerBound, upperBound, digitAccuracy):
    tolerance = 0.5 * 10**(-digitAccuracy)
    midPoint = 0
    error = 0
    fx = lambdify(variable, function)
    iterations = 0
    
    if fx(lowerBound)*fx(upperBound) < 0:
        while ((float)(lowerBound - upperBound) / 2 > tolerance or (float)(lowerBound - upperBound) / 2 < -tolerance) and fx(midPoint) != 0:
            oldMidPoint = midPoint
            midPoint = (float)(lowerBound + upperBound) / 2
            error = (midPoint - oldMidPoint) / midPoint
            if fx(lowerBound) * fx(midPoint) < 0:
                upperBound = midPoint
            else:
                lowerBound = midPoint
            iterations += 1
            
        return iterations, round(midPoint, digitAccuracy), abs(round(error, digitAccuracy))
    else:
        return iterations, None, None
    
def bisectionMethodIntervalLength(lowerBound, upperBound, iterations, digitAccuracy):
    return abs(round((float)(upperBound - lowerBound) / (2**iterations), digitAccuracy))

def bisectionMethodConvergenceRate(variable, function, approximateRoot, digitAccuracy):
    dfx = lambdify(variable, function.diff(variable))
    
    return abs(round(dfx(approximateRoot), digitAccuracy))

def fixedPointIteration (digitAccuracy, guess, variable, functionList):
    tolerance = 0.5 * 10**(-digitAccuracy)
    lastGuess = 0
    x = guess
    usedFunction = None
    
    try:
        for function in functionList:
            dfx = lambdify(variable,function.diff(variable))
            if "j" not in str(dfx(x)) and abs(dfx(x)) < 1:
                if usedFunction != None:
                    oldDFx = lambdify(variable, sympify(usedFunction).diff(variable))
                    
                    if dfx(x) < oldDFx(x):
                        usedFunction = function
                else:
                    usedFunction = function
    except NotImplementedError as e:
        usedFunction = functionList[0]
        
    if usedFunction == None:
        usedFunction = functionList[0]
        
    fx = lambdify(variable, usedFunction)
    iterations = 0
    
    try:
        if lambdify(variable,usedFunction.diff(variable))(x) < 1:
            while abs(x - lastGuess) > tolerance:
                lastGuess = x
                x = fx(lastGuess)
                
                if abs(x - lastGuess) <= tolerance:
                    return iterations, round(x, digitAccuracy), usedFunction
                
                iterations += 1
    except OverflowError as e:
        return -1, None, None
    
    return -1, None, None
    
def newtonsMethod (digitAccuracy, guess, realRoot, variable, function):
    x = guess
    iterations = 1
    error = 0.5 * 10**(-digitAccuracy)
    errorQuadraticSum = 0
    errorLocalSum = 0
    
    while abs(x) > error:
        fx = lambdify(variable, function)
        dfx = lambdify(variable, function.diff(variable))
        
        if forwardError(realRoot, x, digitAccuracy) != 0:
            errorQuadraticSum += float(forwardError(realRoot, x - fx(x)/dfx(x), digitAccuracy)/(forwardError(realRoot, x, digitAccuracy)**2))
            errorLocalSum += float(forwardError(realRoot, x - fx(x)/dfx(x), digitAccuracy)/forwardError(realRoot, x, digitAccuracy))
            
        if abs(fx(x)) < error:
            return iterations, round(x, digitAccuracy), round(errorQuadraticSum, digitAccuracy), round(errorLocalSum, digitAccuracy)
        
        if dfx(x) == 0:
            return iterations, x, 0, 0
        
        x = x - fx(x)/dfx(x)
        iterations += 1
        
    return iterations, x, 0, 0

def forwardError(actualRoot, approximation, digitAccuracy):
    if type(actualRoot) is list:
        smallestDifference = abs(actualRoot[0] - approximation)
        
        for root in actualRoot:
            if abs(root - approximation) < smallestDifference:
                smallestDifference = abs(root - approximation)
                
        return abs(round(smallestDifference, digitAccuracy))
    else:
        return abs(round(abs(actualRoot - approximation), digitAccuracy))

def backwardError(variable, function, approximation, digitAccuracy):
    fx = lambdify(variable, function)
    return abs(round(fx(approximation), digitAccuracy))

functionVariable = Symbol('x')
functionList = []

#functionList = [['x^3 + 2*x + 2', 0.5, -0.770917, -1, 0] , ['E^x + x - 7', 1.5, 1.67282, 1, 2], 
#                ['E^x + sin(x) - 4', 1.5, 1.12998, 1, 2]]


yesList = ["y", "yes"]
noList = ["n", "no"]

digitAccuracy = -1
counter = 0
inputFinished = False

while digitAccuracy < 0:
    digitAccuracy = input("Please enter the number of digits after the decimal for accuracy: ")
    
    if digitAccuracy.isdigit() and int(digitAccuracy) >= 0:
        digitAccuracy = int(digitAccuracy)
    else:
        digitAccuracy = -1
        
while not inputFinished:
    try:
        newFunctionString = input("Please enter a valid function you would like to estimate: ").replace("e", "E")
        newFunction = sympify(newFunctionString)
        lowerBound = float(input("Please enter a lower bound to iterate through: "))
        upperBound = float(input("Please enter an upper bound to iterate through: "))
        guess = float(input("Please enter an initial guess for iteration: "))
        
        try:
            solution = solve(newFunction, functionVariable)
        except NotImplementedError as e:
            solution = nsolve(newFunction, functionVariable, guess)
          
        if type(solution) is list:
            solution = [x for x in solution if "I" not in str(x)]
            if len(solution) == 1:
                solution = solution[0]
            else:
                print()
            
        if lowerBound <= guess <= upperBound and lowerBound != upperBound:
            functionList.append([newFunctionString, guess, solution, lowerBound, upperBound])
        else:
            if lowerBound > upperBound:
                print("Error, lower bound is higher than upper bound.")
            elif lowerBound == upperBound:
                print("Error, lower bound is the same as the upper bound.")
            elif guess < lowerBound or guess > upperBound:
                print("Error, guess is outside the given range")
                
            inputFinished = "Error"
    except SympifyError as e:
        print("Error, invalid input.")
        inputFinished = "Error"
    except ValueError as e:
        print("Error, invalid input.")
        inputFinished = "Error"
    
    if inputFinished != "Error":
        inputFinished = None
        while inputFinished == None:
            inputFinished = input("Would you like to enter more functions to evaluate (Y/N)?")
            
            if inputFinished.lower() in yesList or inputFinished.lower() in noList:
                inputFinished = inputFinished.lower() in noList
            else:
                inputFinished = None
    else:
        inputFinished = False
            
with open("output.txt", "w") as output:
    output.write("Bisection Method:\n")
    while counter < len(functionList):
        result = BisectionMethodFunctionResult(functionList[counter][0], functionList[counter][1],
                                functionList[counter][3], functionList[counter][4])
        
        result.iterations, result.approximateRoot, result.error = bisectionMethod(functionVariable, 
                                          result.function, result.lowerBound, 
                                          result.upperBound, digitAccuracy)
        
        if result.approximateRoot == None:
            output.write("Expression " + str(counter + 1) + ": " + result.functionString)
            output.write("\nCannot approximate root with the given bounds.")
        else:
            result.intervalLength = bisectionMethodIntervalLength(result.lowerBound, 
                                              result.upperBound, result.iterations,
                                              digitAccuracy)
            
            result.convergenceRate = round(abs(1/(2^result.iterations)), digitAccuracy)
            
            output.write("Expression " + str(counter + 1) + ": " + result.functionString)
            output.write("\nApproximate root: " + str(result.approximateRoot))
            output.write("\nIterations needed: " + str(result.iterations))
            output.write("\nInterval Length: " + str(result.intervalLength))
            output.write("\nError: +- " + str(result.error))
            output.write("\nConvergence Type: " + ("Converges" if result.convergenceRate < 1 else "Diverges"))
            output.write("\nConvergence Rate: " + str(result.convergenceRate))
                
        output.write("\n\n")
        counter += 1
    output.write("Fixed Point Iteration:\n")
    
    counter = 0
    
    while counter < len(functionList):
        variationList = []
        useVariations = None
        usedVariation = ""
        
        result = StandardFunctionResult(functionList[counter][0], 
                                functionList[counter][1], functionList[counter][2])
        if result.functionString.count("x") > 1:
            print("\n" + result.functionString)
            for x in range(result.functionString.count("x")):
                usedVariation = input("Please enter a valid function as a variation of the above formula: ")
                
                try:
                    variationList.append(sympify(usedVariation.replace("e", "E")))
                except SympifyError as e:
                    print("Error, invalid input.")
                    x -= 1
        else:
            variationList.append(result.function)
                
        result.iterations, result.approximateRoot, usedVariation = fixedPointIteration(digitAccuracy, result.initialGuess, 
                                          functionVariable, variationList)
        
        output.write("Expression " + str(counter + 1) + ": " + result.functionString)
        output.write("\nUsed variation: " + str(usedVariation))
        
        if result.approximateRoot != None:
            dfx = lambdify(functionVariable, usedVariation.diff(functionVariable))
            
            result.forwardError = forwardError(result.realRoot,
                                              result.approximateRoot, digitAccuracy)
            
            result.backwardsError = backwardError(functionVariable, 
                                              usedVariation, result.approximateRoot,
                                              digitAccuracy)
            result.convergenceRate = round(abs(dfx(result.approximateRoot)), digitAccuracy)
            
            output.write("\nApproximate root: " + str(result.approximateRoot))
            output.write("\nIterations needed: " + str(result.iterations))
            output.write("\nForward Error: " + str(result.forwardError))
            output.write("\nBackward Error: " + str(result.backwardsError))
            output.write("\nConvergence Type: " + ("Converges" if result.convergenceRate < 1 else "Diverges"))
            output.write("\nConvergence Rate: " + str(result.convergenceRate))
        else:
            output.write("\nDoes Not Converge With The Given Function")
            
        output.write("\n\n")
        counter += 1
    output.write("Newton's Method:\n")
    
    counter = 0
    quadraticRate = 0
    localRate = 0
    
    while counter < len(functionList):
        result = StandardFunctionResult(functionList[counter][0], 
                                functionList[counter][1], functionList[counter][2])
        
        result.iterations, result.approximateRoot, quadraticRate, localRate = newtonsMethod(digitAccuracy, 
                                          result.initialGuess, result.realRoot, functionVariable, result.function)
        
        result.forwardError = forwardError(result.realRoot,
                                          result.approximateRoot, digitAccuracy)
        
        result.backwardsError = backwardError(functionVariable, 
                                          result.function, result.approximateRoot,
                                          digitAccuracy)
        
        output.write("Expression " + str(counter + 1) + ": " + result.functionString)
        output.write("\nApproximate root: " + str(result.approximateRoot))
        output.write("\nIterations needed: " + str(result.iterations))
        output.write("\nForward Error: " + str(result.forwardError))
        output.write("\nBackward Error: " + str(result.backwardsError))
        
        fx = lambdify(functionVariable, result.function)
        dfx = lambdify(functionVariable, result.function.diff(functionVariable))
        
        if round(fx(result.approximateRoot), digitAccuracy) == 0 and round(dfx(result.approximateRoot), digitAccuracy) != 0:
            output.write("\nQuadratic Convergence Rate: " + str(quadraticRate))
        else:
            output.write("\nQuadratic Convergence Rate: N/A")
        
        output.write("\nLocal Convergence Rate: " + str(localRate))
            
        output.write("\n\n")
        counter += 1
        
output.close()
