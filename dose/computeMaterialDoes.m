function computeMaterialDoes( map, doseTotal )

% print additional information
doseTissue = doseTotal .* map.mapTissue;
meanDoseTissue = sum(doseTissue(:)) / sum(map.mapTissue(:));
maxDoseTissue = max(doseTissue(:));
doseBone = doseTotal .* map.mapBone;
meanDoseBone = sum(doseBone(:)) / sum(map.mapBone(:));
maxDoseBone = max(doseBone(:));
doseFillings = doseTotal .* map.mapFillings;
meanDoseFillings = sum(doseFillings(:)) / sum(map.mapFillings(:));


% maxDoseFillings = max(doseFillings(:));
fprintf('Mean dose in tissue:       %5.3fGy\n', meanDoseTissue);
fprintf('Mean dose in bone:         %5.3fGy\n', meanDoseBone);
fprintf('Mean dose in fillings:     %5.3fGy\n', meanDoseFillings);
fprintf('Max dose in tissue:        %5.3fGy\n', maxDoseTissue);
fprintf('Max dose in bone:          %5.3fGy\n', maxDoseBone);





end