function[reflection, transmission, absorption] = RTA_iter(energyStep, energy, damping, potentialPositions, potentialValues, delta, energyScale)
    %% Function to calculate reflection, transmission and absorption coefficients
    % This function calculates the reflection, transmission and absorption coefficients
    % for a given energy step, energy, damping, potential positions, potential values,
    % delta and energy scale.
    %
    % Parameters:
    % energyStep: The energy step
    % energy: The energy
    % damping: The damping
    % potentialPositions: The potential positions
    % potentialValues: The potential values
    % delta: The delta
    % energyScale: The energy scale
    %
    % Returns:
    % reflection: The reflection coefficient
    % transmission: The transmission coefficient
    % absorption: The absorption coefficient
        
    % Get the number of potential positions
    [~, numberOfPotentialPositions] = size(potentialPositions); 
    
    % Calculate the initial wave number
    waveNumberInitial = sqrt((energy + 1i * damping) / energyScale); 
    
    % Calculate the final wave number
    waveNumberFinal = sqrt((energy + 1i * damping - energyStep) / energyScale); 
    
    % Calculate the transmission coefficient
    transmissionCoefficient = 2 * waveNumberInitial / (waveNumberInitial + waveNumberFinal); 
    
    % Calculate the reflection coefficient
    reflectionCoefficient = (waveNumberInitial - waveNumberFinal) / (waveNumberFinal + waveNumberInitial); 
    
    % Initialize the potential matrix
    potentialMatrix = zeros(numberOfPotentialPositions, numberOfPotentialPositions); 
    
    for ii = 1:numberOfPotentialPositions
        % Function to the right of the step [A.65]
        waveFunctionRightOfStep(ii) = transmissionCoefficient * exp(1i * waveNumberFinal * potentialPositions(ii)); 
        
        % Equation[4.39]
        potentialMatrix(ii, ii) = delta * potentialValues(ii) / energyScale; 
        
        for jj = 1:numberOfPotentialPositions
            % Equation [4.39]
            GreensFununction(ii, jj) = GreensFun(energyStep, potentialPositions(ii), potentialPositions(jj), energy, damping, energyScale); 
        end
    end
    
    % Calculate the scattering matrix
    scatteringMatrix = eye(numberOfPotentialPositions) - GreensFununction * potentialMatrix; 

    % Solution inside the perturbated region [4.51]
    waveFunctionSolution = scatteringMatrix \ (waveFunctionRightOfStep.'); 
    
    % The first part comes from a pre-existing reflection in equation [A.65]
    reflectedWaveFunction = reflectionCoefficient * exp(-1i * waveNumberInitial * (potentialPositions(1) - 1)) + ExtraTerm(potentialPositions(1) - 1, potentialPositions, potentialValues, delta, waveFunctionSolution, energyStep, energy, damping, energyScale); 
    
    % The first part comes from an already existing transmission in equation [A.65]
    transmittedWaveFunction = transmissionCoefficient * exp(1i * waveNumberFinal * (potentialPositions(numberOfPotentialPositions) + 1)) + ExtraTerm(potentialPositions(numberOfPotentialPositions) + 1, potentialPositions, potentialValues, delta, waveFunctionSolution, energyStep, energy, damping, energyScale); 
    
    % Calculate the transmission
    transmission = (real(waveNumberFinal) / real(waveNumberInitial)) * abs(transmittedWaveFunction)^2; 
    
    % Compare the reflected part with just an incoming wave
    reflection = (abs(reflectedWaveFunction / (exp(1i * waveNumberInitial * (potentialPositions(1) - 1)))))^2; 
    
    % Calculate the absorption
    absorption = 1 - (reflection + transmission); 
end