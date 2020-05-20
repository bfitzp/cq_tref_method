function b_lr = InterpolateIncField(z_lr,Z_hr,b_hr)

diff_lr = abs(z_lr - z_lr(1,:));
diff_hr = abs(Z_hr - Z_hr(1,:));

b_lr = interp1(diff_hr,b_hr,diff_lr);

end