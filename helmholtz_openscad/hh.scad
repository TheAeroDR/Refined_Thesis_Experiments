// hh.scad
// Derive coil radii from required uniform-field radius

// ========== PARAMETERS ==========
uniform_radius = 40;
channel_width = 6;
channel_height = 4;
channel_base = 1;
strut_width = channel_height-channel_base;
strut_depth = 4;
$fn = 200;
clearance = 2;

coil_radius_1 = uniform_radius / 0.31;
coil_radius_2 = coil_radius_1 + channel_width*2 + clearance;
coil_radius_3 = coil_radius_2 + channel_width*2 + clearance;

sep_1 = coil_radius_1;
sep_2 = coil_radius_2;
sep_3 = coil_radius_3;

// ========== MODULES ==========
module U_channel(r, width=channel_width, height=channel_height, base=channel_base) {
    x_offset = r-base;
    y_offset = -width/2;
    pts = [
        [x_offset, y_offset],
        [height + x_offset, y_offset],
        [height + x_offset, base + y_offset],
        [base + x_offset, base + y_offset],
        [base + x_offset, width - base + y_offset],
        [height + x_offset, width - base + y_offset],
        [height + x_offset, width + y_offset],
        [x_offset, width + y_offset]
    ];
    rotate_extrude($fn=$fn) polygon(pts);
}

module HH_pair(r, separation) {
translate([0,0, separation/2]) U_channel(r);
translate([0,0,-separation/2]) U_channel(r);
}
module pair_struts(r, separation) {
for(angle=[0:60:359]) {
rotate([0,0,angle])
translate([r,0,-separation/2])
cube([strut_width, strut_depth, separation]);
}
}
module axis_X(r, separation) {
rotate([0,90,0]) {
HH_pair(r, separation);
pair_struts(r, separation-channel_width);
}
}
module axis_Y(r, separation) {
rotate([90,0,0]) {
HH_pair(r, separation);
pair_struts(r, separation-channel_width);
}
}
module axis_Z(r, separation) {
HH_pair(r, separation);
pair_struts(r, separation-channel_width);
}
// ========== RENDER =============
axis_X(coil_radius_1, sep_1);
axis_Y(coil_radius_2, sep_2);
axis_Z(coil_radius_3, sep_3);
