function[extraTermSolution] = ExtraTerm(position, potentialPositions, potentialValues, delta, waveFunctionSolution, energyStep, energy, damping, energyScale)
    %% Function to calculate the extra term in the wave function solution
    % This function calculates the extra term in the wave function solution
    % for a given position, potential positions, potential values, delta, wave function solution,
    % energy step, energy, damping, and energy scale. It does this by iterating over each potential position
    % and summing the product of the Green's function, scale, potential value, and wave function solution.
    %
    % Parameters:
    % position: The position
    % potentialPositions: The potential positions
    % potentialValues: The potential values
    % delta: The delta
    % waveFunctionSolution: The wave function solution
    % energyStep: The energy step
    % energy: The energy
    % damping: The damping
    % energyScale: The energy scale
    %
    % Returns:
    % extraTermSolution: The extra term in the wave function solution
    
        extraTermSolution = 0; % Initialize the extra term solution
        [~, numPotentialPositions] = size(potentialPositions); % Get the number of potential positions
        scale = delta / energyScale; % Calculate the scale
    
        % Iterate over each potential position
        parfor i = 1:numPotentialPositions
            % Calculate the extra term in the wave function solution (equation [4.36])
            extraTermSolution = extraTermSolution + GreensFun(energyStep, position, potentialPositions(i), energy, damping, energyScale) * scale * potentialValues(i) * waveFunctionSolution(i);
        end
    end