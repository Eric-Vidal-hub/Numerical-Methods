%% OUTSIDE SOLUTION
% Function to calculate the extra part added to the solution due to the potential.
%
% Syntax:  extraSolution = ExtraTerm(position, potentialPositions, potentialValues, delta, waveFunctionSolution, energyStep, energy, damping, energyScale)
%
% Inputs:
%    position - The position for which to calculate the extra solution
%    potentialPositions - The positions of the potential
%    potentialValues - The values of the potential at the corresponding positions
%    delta - A small change in position
%    waveFunctionSolution - The solution of the wave function at the potential positions
%    energyStep - A parameter related to the energy level
%    energy - The energy level
%    damping - A damping factor
%    energyScale - A scaling factor for the energy
%
% Outputs:
%    extraSolution - Calculated extra part of the solution

function[extraSolution]= ExtraTerm(position, potentialPositions, potentialValues, delta, waveFunctionSolution, energyStep, energy, damping, energyScale)
    extraSolution=0;
    [~,numberOfPotentialPositions]=size(potentialPositions);
    scale = delta/energyScale;
    for ii=1:numberOfPotentialPositions
        extraSolution=extraSolution+GreensFun(energyStep, position, potentialPositions(ii), energy, damping, energyScale)*scale*potentialValues(ii)*waveFunctionSolution(ii); %equation [4.36]
    end
end