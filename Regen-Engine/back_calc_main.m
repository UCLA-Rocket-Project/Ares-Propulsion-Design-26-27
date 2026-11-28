clc; clear; close all;

C = 1;
known_temp_station = 100;
known_temp = 250;
Knockdown.is_found = false;
temp_tol = 1;
temp_error = realmax;

while (abs(temp_error) > temp_tol)
    Knockdown = regen_back_calc(C);
    calc_temp = Arrays.T_bulk_array(known_temp_station);
    temp_error = calc_temp - known_temp;
    temp_range = abs(Arrays.T_bulk_array(1) - Arrays.T_bulk_array(end));
    if (abs(temp_error) > temp_tol)
        C = C + 0.3*(temp_error/temp_range);
    end
end

print(C);