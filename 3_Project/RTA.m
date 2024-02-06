function[reflection, transmission, absorption] = RTA(energyStep, energy, damping, potentialPositions, potentialValues, delta, energyScale)
    %% Function to calculate reflection, transmission and absorption coefficients
    % This function calculates the reflection, transmission and absorption coefficients
    % for a given energy step, energy, damping, potential positions, potential values,
    % delta and energy scale. It does this by iterating over each energy value and
    % calculating the reflection, transmission, and absorption coefficients using the
    % RTA_iter function. It then plots these coefficients against the energy values.
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
        numE=length(energy); % Get the number of energies
        for kk=1:numE
            [reflection(kk), transmission(kk), absorption(kk)] = RTA_iter(energyStep, energy(kk), damping, potentialPositions, potentialValues, delta, energyScale); % Calculate the reflection and transmission
        end
    
        % Plot of R(E), T(E) and A(E)
        figure;
        plot(energy, reflection, 'LineWidth', 2);
        hold on;
        plot(energy, transmission, 'LineWidth', 2);
        plot(energy, absorption, 'LineWidth', 2);
        xlabel('E (eV)','FontSize',26);
        ylabel('R(E), T(E), A(E)','FontSize',26);
        fontsize(gca, 22,'points')   % 'pixels', 'centimeters', 'inches'
        ylim([0,1]);
        set(gca,'Box','on');
        set(gca,'linewidth',1);
        grid on;
        % Set the color of the axes and the grid to a light gray
        set(gca, 'Color', [0.9 0.9 0.9], 'GridColor', [0.5 0.5 0.5]);
        % Set the background color of the figure to a darker gray
        set(gcf, 'Color', [0.7 0.7 0.7]);
        % set(gca,'TickLength',[0.03, 0.02]);
        legend('R(E)', 'T(E)', 'A(E)', 'Location', 'Best');
    end