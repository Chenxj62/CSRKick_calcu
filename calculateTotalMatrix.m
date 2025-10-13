function R_total = calculateTotalMatrix(beamline)
    % 计算整个束线的一阶传输矩阵
    R_total = eye(6);
    for i = 1:length(beamline)
        R_elem = getElementMatrix(beamline(i));
        R_total = R_elem * R_total;
    end
end

