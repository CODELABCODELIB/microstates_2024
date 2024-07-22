function [boot_data] = create_boot_test_val(test_val)
n_boots = 1000;
boot_data = zeros(2500, n_boots);
for boot=1:n_boots
    boot_data(:,boot) = datasample(test_val,length(test_val));
end
end