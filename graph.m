% Animates the snapshots written by kinetic-solver (MATLAB or GNU Octave).
%
% Usage, from the repository root:
%     graph            % reads the default "output" directory
%     graph('runs/dam-break')
%
% Each snapshot_NNNNN.dat file has the columns x z h hu [h_exact hu_exact];
% times.dat lists "index t" for every snapshot.

function graph(directory)
	if nargin < 1
		directory = 'output';
	end

	times = load(fullfile(directory, 'times.dat'));
	for k = 1:size(times, 1)
		file = fullfile(directory, sprintf('snapshot_%05d.dat', times(k, 1)));
		data = dlmread(file, ' ', 2, 0);   % skip the two '#' header lines
		x = data(:, 1);
		z = data(:, 2);
		h = data(:, 3);
		hu = data(:, 4);

		subplot(2, 1, 1);
		area_bottom = min(z) - 0.05;
		plot(x, h + z, 'b', x, z, 'k');
		hold on
		if size(data, 2) >= 6
			plot(x, data(:, 5) + z, 'r--');
		end
		hold off
		ylim([area_bottom, max(h + z) + 0.1]);
		title(sprintf('t = %.3f s', times(k, 2)));
		ylabel('elevation (m)');

		subplot(2, 1, 2);
		plot(x, hu, 'b');
		if size(data, 2) >= 6
			hold on
			plot(x, data(:, 6), 'r--');
			hold off
		end
		ylabel('hu (m^2/s)');
		xlabel('x (m)');
		drawnow;
		pause(0.05);
	end
end
