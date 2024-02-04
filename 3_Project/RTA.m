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
function[reflection, transmission, absorption] = RTA(energyStep, energy, damping, potentialPositions, potentialValues, delta, energyScale)
    [~, numberOfPotentialPositions] = size(potentialPositions); % Get the number of potential positions
    waveNumberInitial = sqrt((energy + 1i * damping) / energyScale); % Calculate the initial wave number
    waveNumberFinal = sqrt((energy + 1i * damping - energyStep) / energyScale); % Calculate the final wave number
    transmissionCoefficient = 2 * waveNumberInitial / (waveNumberInitial + waveNumberFinal); % Calculate the transmission coefficient
    reflectionCoefficient = (waveNumberInitial - waveNumberFinal) / (waveNumberFinal + waveNumberInitial); % Calculate the reflection coefficient
    potentialMatrix = zeros(numberOfPotentialPositions, numberOfPotentialPositions); % Initialize the potential matrix
    for ii = 1:numberOfPotentialPositions
     waveFunctionRightOfStep(ii) = transmissionCoefficient * exp(1i * waveNumberFinal * potentialPositions(ii)); % Function to the right of the step [A.65]
     potentialMatrix(ii, ii) = delta * potentialValues(ii) / energyScale; % Equation[4.39]
     for jj = 1:numberOfPotentialPositions
         GreensFununction(ii, jj) = GreensFun(energyStep, potentialPositions(ii), potentialPositions(jj), energy, damping, energyScale); % Equation [4.39]
     end
    end
    scatteringMatrix = eye(numberOfPotentialPositions) - GreensFununction * potentialMatrix; % Calculate the scattering matrix
 
    waveFunctionSolution = scatteringMatrix \ (waveFunctionRightOfStep.'); % Solution inside the perturbated region [4.51]
    reflectedWaveFunction = reflectionCoefficient * exp(-1i * waveNumberInitial * (potentialPositions(1) - 1)) + ExtraTerm(potentialPositions(1) - 1, potentialPositions, potentialValues, delta, waveFunctionSolution, energyStep, energy, damping, energyScale); % The first part comes from a pre-existing reflection in equation [A.65]
    transmittedWaveFunction = transmissionCoefficient * exp(1i * waveNumberFinal * (potentialPositions(numberOfPotentialPositions) + 1)) + ExtraTerm(potentialPositions(numberOfPotentialPositions) + 1, potentialPositions, potentialValues, delta, waveFunctionSolution, energyStep, energy, damping, energyScale); % The first part comes from an already existing transmission in equation [A.65]
    
    transmission = (real(waveNumberFinal) / real(waveNumberInitial)) * abs(transmittedWaveFunction)^2; % Calculate the transmission
    reflection = (abs(reflectedWaveFunction / (exp(1i * waveNumberInitial * (potentialPositions(1) - 1)))))^2; % Compare the reflected part with just an incoming wave
    absorption = 1 - (reflection + transmission); % Calculate the absorption
 end