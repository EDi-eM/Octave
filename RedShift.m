function redshift_simulator()
    % ANIMATED REDSHIFT EFFECT SIMULATOR - Octave Compatible
    % Simulates Doppler, Cosmological, and Gravitational redshift effects with animations

    clear; clc; close all;

    % Check if running on Octave
    if exist('OCTAVE_VERSION', 'builtin')
        pkg load signal;  % Load signal package if available
        graphics_toolkit('qt'); % Use Qt graphics toolkit for better performance
    end

    % Physical constants
    c = 299792458; % Speed of light (m/s)
    H0 = 70; % Hubble constant (km/s/Mpc)

    % Create main figure
    figure('Position', [100, 100, 1200, 800], 'Name', 'Animated Redshift Effect Simulator');

    % Menu for different redshift types
    fprintf('=== ANIMATED REDSHIFT EFFECT SIMULATOR ===\n');
    fprintf('1. Doppler Redshift Animation (Moving Object)\n');
    fprintf('2. Cosmological Redshift Animation (Expanding Universe)\n');
    fprintf('3. Gravitational Redshift Animation (Orbiting Light)\n');
    fprintf('4. Spectral Line Animation (All Effects)\n');
    fprintf('5. Simple Interactive Demo\n');

    choice = input('Select simulation type (1-5): ');

    switch choice
        case 1
            animate_doppler_redshift();
        case 2
            animate_cosmological_redshift();
        case 3
            animate_gravitational_redshift();
        case 4
            animate_spectral_lines();
        case 5
            simple_interactive_demo();
        otherwise
            fprintf('Invalid choice. Running Doppler animation...\n');
            animate_doppler_redshift();
    end
end

function animate_doppler_redshift()
    % Animate Doppler redshift with moving object

    clf;
    c = 299792458; % m/s
    lambda0 = 656.3; % Hydrogen alpha line (nm)

    % Animation parameters
    t_max = 10; % seconds
    dt = 0.2; % Slower for Octave
    t_steps = 0:dt:t_max;

    % Object motion parameters
    v_max = 0.8 * c; % Maximum velocity (80% speed of light)

    fprintf('Starting Doppler Redshift Animation...\n');
    fprintf('Object oscillates between -0.8c and +0.8c\n');
    fprintf('Red = moving away, Blue = moving toward observer\n');
    fprintf('Press Ctrl+C to stop animation\n');

    for i = 1:length(t_steps)
        t = t_steps(i);
        clf;

        % Current velocity (sinusoidal motion)
        v_current = v_max * sin(2*pi*t/t_max);
        v_fraction = v_current / c;

        % Calculate relativistic Doppler shift
        gamma = 1 / sqrt(1 - v_fraction^2);
        lambda_observed = lambda0 * gamma * (1 + v_fraction) / sqrt(1 - v_fraction^2);
        z_current = (lambda_observed - lambda0) / lambda0;

        % Color based on redshift - simplified for Octave
        if z_current > 0
            obj_color = 'r'; % Red for redshift
        else
            obj_color = 'b'; % Blue for blueshift
        end

        % Main animation plot
        subplot(2,3,1);

        % Draw observer
        plot(0, 0, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', 'black');
        hold on;

        % Draw moving object
        x_pos = 5 + 2*sin(2*pi*t/t_max); % Oscillating position
        plot(x_pos, 0, 'o', 'MarkerSize', 20, 'Color', obj_color, 'MarkerFaceColor', obj_color);

        % Draw light wave - simplified
        wave_x = linspace(0, x_pos, 50);
        if x_pos > 0 % Object to the right
            wave_freq = 1 / (lambda_observed/lambda0); % Frequency scaling
            wave_y = 0.3 * sin(2*pi*wave_freq*(wave_x - c*t/1e8));
            plot(wave_x, wave_y, obj_color, 'LineWidth', 2);
        end

        % Velocity arrow - using plot instead of quiver for compatibility
        arrow_scale = v_fraction * 1.5;
        if arrow_scale > 0
            plot([x_pos, x_pos + arrow_scale], [-0.5, -0.5], 'r-', 'LineWidth', 3);
            plot(x_pos + arrow_scale, -0.5, 'r>', 'MarkerSize', 8);
        elseif arrow_scale < 0
            plot([x_pos, x_pos + arrow_scale], [-0.5, -0.5], 'r-', 'LineWidth', 3);
            plot(x_pos + arrow_scale, -0.5, 'r<', 'MarkerSize', 8);
        end

        xlim([-1, 8]);
        ylim([-1, 1]);
        xlabel('Distance');
        title(sprintf('Doppler Effect Animation - t = %.1fs', t));
        text(1, 0.7, sprintf('v = %.2fc', v_fraction), 'FontSize', 12);
        text(1, 0.5, sprintf('lambda = %.1f nm', lambda_observed), 'FontSize', 12);
        text(1, 0.3, sprintf('z = %.3f', z_current), 'FontSize', 12);
        grid on;

        % Velocity vs time plot
        subplot(2,3,2);
        plot(t_steps(1:i), v_max*sin(2*pi*t_steps(1:i)/t_max)/c, 'b-', 'LineWidth', 2);
        hold on;
        plot(t, v_fraction, 'ro', 'MarkerSize', 8);
        xlabel('Time (s)');
        ylabel('Velocity (fraction of c)');
        title('Velocity vs Time');
        xlim([0, t_max]);
        ylim([-1, 1]);
        grid on;

        % Redshift vs time plot
        subplot(2,3,3);
        z_history = zeros(1, i);
        for j = 1:i
            v_temp = v_max*sin(2*pi*t_steps(j)/t_max)/c;
            lambda_temp = lambda0 * (1 + v_temp) / sqrt(1 - v_temp^2);
            z_history(j) = (lambda_temp - lambda0) / lambda0;
        end
        plot(t_steps(1:i), z_history, 'g-', 'LineWidth', 2);
        hold on;
        plot(t, z_current, 'ro', 'MarkerSize', 8);
        xlabel('Time (s)');
        ylabel('Redshift z');
        title('Redshift vs Time');
        xlim([0, t_max]);
        grid on;

        % Spectrum plot
        subplot(2,3,[4,5]);
        frequencies = linspace(400, 800, 200); % Reduced points for speed
        spectrum_rest = exp(-((frequencies - lambda0)/10).^2);
        spectrum_observed = exp(-((frequencies - lambda_observed)/10).^2);

        plot(frequencies, spectrum_rest, 'k--', 'LineWidth', 1.5);
        hold on;
        plot(frequencies, spectrum_observed, obj_color, 'LineWidth', 2);
        xlabel('Wavelength (nm)');
        ylabel('Intensity');
        title('Spectral Line Shift');
        legend('Rest wavelength', 'Observed', 'Location', 'northeast');
        xlim([400, 800]);
        grid on;

        % Progress bar
        subplot(2,3,6);
        progress = i / length(t_steps);
        barh(1, progress, 'FaceColor', 'green');
        xlim([0, 1]);
        ylim([0.5, 1.5]);
        xlabel('Animation Progress');
        title(sprintf('%.1f%% Complete', progress*100));

        drawnow;
        pause(dt);
    end

    fprintf('Animation complete!\n');
end

function animate_cosmological_redshift()
    % Animate cosmological redshift in expanding universe

    clf;
    H0 = 70; % km/s/Mpc
    c_km = 299792.458; % km/s
    lambda0 = 656.3; % nm

    % Animation parameters
    t_max = 13.8; % Billion years (age of universe)
    dt = 0.5; % Slower for Octave
    t_steps = 0.1:dt:t_max;

    % Galaxy positions and distances
    n_galaxies = 6; % Reduced for performance
    galaxy_distances = linspace(50, 400, n_galaxies); % Mpc
    galaxy_angles = linspace(0, 2*pi, n_galaxies+1);
    galaxy_angles = galaxy_angles(1:end-1);

    fprintf('Starting Cosmological Redshift Animation...\n');
    fprintf('Simulating universe expansion over 13.8 billion years\n');
    fprintf('Galaxies move away - redshift increases with distance\n');
    fprintf('Press Ctrl+C to stop animation\n');

    for i = 1:length(t_steps)
        t = t_steps(i);
        clf;

        % Scale factor (simplified)
        scale_factor = t / t_max;
        current_distances = galaxy_distances * scale_factor;

        % Main universe plot
        subplot(2,2,1);

        % Draw observer at center
        plot(0, 0, 'ko', 'MarkerSize', 15, 'MarkerFaceColor', 'yellow');
        hold on;

        % Draw galaxies with simplified colors
        colors = {'r', 'g', 'b', 'm', 'c', 'k'};
        for g = 1:n_galaxies
            x_pos = current_distances(g) * cos(galaxy_angles(g));
            y_pos = current_distances(g) * sin(galaxy_angles(g));

            % Calculate redshift
            recession_velocity = H0 * current_distances(g);
            z_cosmo = recession_velocity / c_km;

            plot(x_pos, y_pos, 's', 'MarkerSize', 12, 'Color', colors{g}, 'MarkerFaceColor', colors{g});

            % Draw velocity vector - simplified
            v_scale = recession_velocity / 2000; % Scale for visualization
            plot([x_pos, x_pos + v_scale*cos(galaxy_angles(g))], ...
                 [y_pos, y_pos + v_scale*sin(galaxy_angles(g))], 'r-', 'LineWidth', 2);
        end

        max_dist = max(current_distances) * 1.2;
        xlim([-max_dist, max_dist]);
        ylim([-max_dist, max_dist]);
        xlabel('Distance (Mpc)');
        ylabel('Distance (Mpc)');
        title(sprintf('Expanding Universe - Age: %.1f Gyr', t));
        axis equal;
        grid on;

        % Hubble diagram
        subplot(2,2,2);
        recession_velocities = H0 * current_distances;
        plot(current_distances, recession_velocities/1000, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
        xlabel('Distance (Mpc)');
        ylabel('Recession Velocity (1000 km/s)');
        title('Hubble Law');
        grid on;

        % Scale factor vs time
        subplot(2,2,3);
        plot(t_steps(1:i), t_steps(1:i)/t_max, 'g-', 'LineWidth', 2);
        hold on;
        plot(t, scale_factor, 'ro', 'MarkerSize', 8);
        xlabel('Time (Gyr)');
        ylabel('Scale Factor');
        title('Universe Expansion');
        xlim([0, t_max]);
        ylim([0, 1.1]);
        grid on;

        % Redshift evolution - simplified
        subplot(2,2,4);
        z_current = H0 * galaxy_distances * scale_factor / c_km;
        plot(current_distances, z_current, 'ro-', 'LineWidth', 2);
        xlabel('Current Distance (Mpc)');
        ylabel('Redshift z');
        title('Distance vs Redshift');
        grid on;

        drawnow;
        pause(dt);
    end

    fprintf('Animation complete!\n');
end

function animate_gravitational_redshift()
    % Animate gravitational redshift for orbiting light source

    clf;
    G = 6.674e-11;
    c = 299792458;
    M_star = 2 * 1.989e30; % 2 solar masses
    lambda0 = 656.3; % nm

    % Orbital parameters
    r_orbit = 50000e3; % 50,000 km orbit radius
    period = 2 * pi * sqrt(r_orbit^3 / (G * M_star)); % Orbital period

    % Animation parameters
    t_max = period;
    dt = period / 50; % Reduced steps for Octave
    t_steps = 0:dt:t_max;

    fprintf('Starting Gravitational Redshift Animation...\n');
    fprintf('Light source orbiting a massive star (2 solar masses)\n');
    fprintf('Redshift varies with distance from gravitational source\n');
    fprintf('Press Ctrl+C to stop animation\n');

    for i = 1:length(t_steps)
        t = t_steps(i);
        clf;

        % Orbital position
        theta = 2 * pi * t / period;
        x_pos = r_orbit * cos(theta);
        y_pos = r_orbit * sin(theta);

        % Distance from star center
        r_current = sqrt(x_pos^2 + y_pos^2);

        % Gravitational redshift
        z_grav = G * M_star / (c^2 * r_current);
        lambda_observed = lambda0 * (1 + z_grav);

        % Main orbital plot
        subplot(2,2,1);

        % Draw central star
        plot(0, 0, 'yo', 'MarkerSize', 25, 'MarkerFaceColor', 'yellow');
        hold on;

        % Draw orbit path
        orbit_theta = linspace(0, 2*pi, 100);
        orbit_x = r_orbit * cos(orbit_theta);
        orbit_y = r_orbit * sin(orbit_theta);
        plot(orbit_x/1000, orbit_y/1000, 'k--', 'LineWidth', 1);

        % Draw orbiting light source
        plot(x_pos/1000, y_pos/1000, 'ro', 'MarkerSize', 15, 'MarkerFaceColor', 'red');

        % Draw light ray to observer
        observer_x = r_orbit * 2;
        observer_y = 0;
        plot(observer_x/1000, observer_y/1000, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'black');

        % Light ray
        plot([x_pos/1000, observer_x/1000], [y_pos/1000, observer_y/1000], 'r-', 'LineWidth', 2);

        xlim([-80, 120]);
        ylim([-80, 80]);
        xlabel('Distance (1000 km)');
        ylabel('Distance (1000 km)');
        title(sprintf('Gravitational Redshift - Phase: %.1f deg', theta*180/pi));
        axis equal;
        grid on;

        % Redshift vs time
        subplot(2,2,2);
        z_history = G * M_star ./ (c^2 * r_orbit); % Constant for circular orbit
        plot(t_steps(1:i)/3600, z_history*1e6*ones(1,i), 'g-', 'LineWidth', 2);
        hold on;
        plot(t/3600, z_grav*1e6, 'ro', 'MarkerSize', 8);
        xlabel('Time (hours)');
        ylabel('Redshift z (x10^-6)');
        title('Gravitational Redshift');
        xlim([0, t_max/3600]);
        grid on;

        % Wavelength shift
        subplot(2,2,3);
        frequencies = linspace(656, 657, 100); % Narrow range for visibility
        spectrum_rest = exp(-((frequencies - lambda0)/0.1).^2);
        spectrum_observed = exp(-((frequencies - lambda_observed)/0.1).^2);

        plot(frequencies, spectrum_rest, 'k--', 'LineWidth', 1.5);
        hold on;
        plot(frequencies, spectrum_observed, 'r-', 'LineWidth', 2);
        xlabel('Wavelength (nm)');
        ylabel('Intensity');
        title('Spectral Line Shift (Magnified)');
        legend('Rest', 'Observed', 'Location', 'northeast');
        grid on;

        % Orbital phase
        subplot(2,2,4);
        phase_degrees = theta * 180 / pi;
        plot(phase_degrees, 1, 'ro', 'MarkerSize', 15);
        hold on;
        plot([0, 360], [1, 1], 'k-', 'LineWidth', 2);
        xlim([0, 360]);
        ylim([0.5, 1.5]);
        xlabel('Orbital Phase (degrees)');
        title(sprintf('Phase: %.1f deg', phase_degrees));

        % Add text information
        text(0.02, 0.95, sprintf('z = %.2e', z_grav), 'Units', 'normalized', 'FontSize', 12);
        text(0.02, 0.90, sprintf('Delta-lambda = %.6f nm', lambda_observed - lambda0), 'Units', 'normalized', 'FontSize', 10);

        drawnow;
        pause(dt);
    end

    fprintf('Animation complete!\n');
end

function animate_spectral_lines()
    % Animate spectral line shifts for all three effects simultaneously

    clf;
    lambda0 = 656.3; % Hydrogen alpha
    c = 299792458;

    % Animation parameters
    t_max = 8;
    dt = 0.2; % Slower for Octave
    t_steps = 0:dt:t_max;

    fprintf('Starting Spectral Line Animation...\n');
    fprintf('Comparing all three redshift effects simultaneously\n');
    fprintf('Press Ctrl+C to stop animation\n');

    for i = 1:length(t_steps)
        t = t_steps(i);
        clf;

        % Time-varying parameters
        v_doppler = 0.6 * c * sin(2*pi*t/t_max); % Oscillating velocity
        d_cosmo = 200 + 100*sin(2*pi*t/(t_max*2)); % Varying distance (Mpc)
        r_grav = 10000e3 + 5000e3*sin(2*pi*t/(t_max*3)); % Varying distance from massive object

        % Calculate redshifts
        z_doppler = (v_doppler/c) / sqrt(1 - (v_doppler/c)^2);
        z_cosmo = 70 * d_cosmo / 299792.458;
        z_grav = 6.674e-11 * 2*1.989e30 / (c^2 * r_grav);

        % Observed wavelengths
        lambda_doppler = lambda0 * (1 + z_doppler);
        lambda_cosmo = lambda0 * (1 + z_cosmo);
        lambda_grav = lambda0 * (1 + z_grav);

        % Spectrum plot
        subplot(2,1,1);
        frequencies = linspace(400, 900, 500); % Reduced for speed

        % Original spectrum
        spectrum_rest = exp(-((frequencies - lambda0)/8).^2);
        plot(frequencies, spectrum_rest, 'k-', 'LineWidth', 3);
        hold on;

        % Shifted spectra
        spectrum_doppler = exp(-((frequencies - lambda_doppler)/8).^2);
        spectrum_cosmo = exp(-((frequencies - lambda_cosmo)/8).^2);
        spectrum_grav = exp(-((frequencies - lambda_grav)/8).^2);

        plot(frequencies, spectrum_doppler + 0.3, 'r-', 'LineWidth', 2);
        plot(frequencies, spectrum_cosmo + 0.6, 'b-', 'LineWidth', 2);
        plot(frequencies, spectrum_grav + 0.9, 'g-', 'LineWidth', 2);

        % Mark line positions
        plot([lambda0, lambda0], [0, 1.2], 'k--', 'LineWidth', 1);
        plot([lambda_doppler, lambda_doppler], [0.3, 1.5], 'r--', 'LineWidth', 1);
        plot([lambda_cosmo, lambda_cosmo], [0.6, 1.8], 'b--', 'LineWidth', 1);
        plot([lambda_grav, lambda_grav], [0.9, 2.1], 'g--', 'LineWidth', 1);

        xlabel('Wavelength (nm)');
        ylabel('Intensity (offset for clarity)');
        title(sprintf('Spectral Line Shifts Animation - t = %.1fs', t));
        legend('Rest (H-alpha)', 'Doppler', 'Cosmological', 'Gravitational', 'Location', 'northeast');
        xlim([620, 720]);
        ylim([0, 2.2]);
        grid on;

        % Information panel
        text(625, 2.0, sprintf('Doppler: v = %.2fc, z = %.4f', v_doppler/c, z_doppler), 'FontSize', 10, 'Color', 'red');
        text(625, 1.8, sprintf('Cosmological: d = %.0f Mpc, z = %.4f', d_cosmo, z_cosmo), 'FontSize', 10, 'Color', 'blue');
        text(625, 1.6, sprintf('Gravitational: r = %.0f km, z = %.2e', r_grav/1000, z_grav), 'FontSize', 10, 'Color', 'green');

        % Redshift comparison plot
        subplot(2,1,2);

        % Current redshift values
        plot(1, z_doppler, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'red');
        hold on;
        plot(2, z_cosmo, 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'blue');
        plot(3, z_grav*1e6, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'green');

        set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'Doppler', 'Cosmological', 'Grav (x10^6)'});
        ylabel('Redshift z');
        title('Current Redshift Values');
        xlim([0.5, 3.5]);
        grid on;

        drawnow;
        pause(dt);
    end

    fprintf('Animation complete!\n');
end

function simple_interactive_demo()
    % Simplified interactive demo for Octave

    clf;
    lambda0 = 656.3;

    fprintf('Simple Interactive Demo\n');
    fprintf('Enter parameters manually:\n');

    while true
        fprintf('\n--- Redshift Calculator ---\n');
        fprintf('1. Doppler effect\n');
        fprintf('2. Cosmological redshift\n');
        fprintf('3. Gravitational redshift\n');
        fprintf('4. Exit\n');

        choice = input('Select option (1-4): ');

        if choice == 4
            break;
        end

        clf;

        switch choice
            case 1
                v_frac = input('Enter velocity as fraction of c (-0.9 to 0.9): ');
                v_frac = max(-0.9, min(0.9, v_frac));
                z = v_frac / sqrt(1 - v_frac^2);
                lambda_obs = lambda0 * (1 + z);

                fprintf('Doppler redshift: z = %.4f\n', z);
                fprintf('Observed wavelength: %.2f nm\n', lambda_obs);

                % Plot spectrum
                freq = linspace(600, 700, 200);
                spec_rest = exp(-((freq - lambda0)/5).^2);
                spec_obs = exp(-((freq - lambda_obs)/5).^2);

                plot(freq, spec_rest, 'k-', 'LineWidth', 2);
                hold on;
                plot(freq, spec_obs, 'r-', 'LineWidth', 2);
                xlabel('Wavelength (nm)');
                ylabel('Intensity');
                title(sprintf('Doppler Shift: v = %.2fc, z = %.4f', v_frac, z));
                legend('Rest', 'Observed', 'Location', 'best');
                grid on;

            case 2
                distance = input('Enter distance in Mpc (10-1000): ');
                distance = max(10, min(1000, distance));
                z = 70 * distance / 299792.458;
                lambda_obs = lambda0 * (1 + z);

                fprintf('Cosmological redshift: z = %.4f\n', z);
                fprintf('Recession velocity: %.0f km/s\n', 70*distance);

                % Plot spectrum
                freq = linspace(600, 800, 200);
                spec_rest = exp(-((freq - lambda0)/5).^2);
                spec_obs = exp(-((freq - lambda_obs)/5).^2);

                plot(freq, spec_rest, 'k-', 'LineWidth', 2);
                hold on;
                plot(freq, spec_obs, 'b-', 'LineWidth', 2);
                xlabel('Wavelength (nm)');
                ylabel('Intensity');
                title(sprintf('Cosmological: d = %.0f Mpc, z = %.4f', distance, z));
                legend('Rest', 'Observed', 'Location', 'best');
                grid on;

            case 3
                mass_solar = input('Enter stellar mass in solar masses (0.5-5): ');
                radius_km = input('Enter distance from star in km (1000-100000): ');

                mass_solar = max(0.5, min(5, mass_solar));
                radius_km = max(1000, min(100000, radius_km));

                G = 6.674e-11;
                M_sun = 1.989e30;
                c = 299792458;

                z = G * mass_solar * M_sun / (c^2 * radius_km * 1000);
                lambda_obs = lambda0 * (1 + z);

                fprintf('Gravitational redshift: z = %.2e\n', z);
                fprintf('Wavelength shift: %.6f nm\n', lambda_obs - lambda0);

                % Plot spectrum (zoomed in)
                freq = linspace(656, 657, 200);
                spec_rest = exp(-((freq - lambda0)/0.05).^2);
                spec_obs = exp(-((freq - lambda_obs)/0.05).^2);

                plot(freq, spec_rest, 'k-', 'LineWidth', 2);
                hold on;
                plot(freq, spec_obs, 'g-', 'LineWidth', 2);
                xlabel('Wavelength (nm)');
                ylabel('Intensity');
                title(sprintf('Gravitational: M = %.1f M_sun, r = %.0f km', mass_solar, radius_km));
                legend('Rest', 'Observed', 'Location', 'best');
                grid on;
        end

        drawnow;
        input('Press Enter to continue...');
    end

    fprintf('Demo ended.\n');
end

% Run the simulator
redshift_simulator();
