/* The "scale" of the system; formally SPACING. */
float ZIPPER_RESOLUTION = 0.0005;

/* mesh display level */
int mesh_level = 3;


init_resolution_parameters()
{
    set_max_edge_length_factor(4.0);

    set_fill_edge_length_factor(2.0);

    set_conf_edge_count_factor(1.0);
    set_conf_edge_zero(0);
    set_conf_angle(0);
    set_conf_exponent(1.0);

    set_align_near_dist_factor(2.0);
    set_align_near_cos(0.3);

    set_eat_near_dist_factor(2.0);
    set_eat_near_cos(-0.5);
    set_eat_start_iters(2);
    set_eat_start_factor(4.0);    

    set_clip_near_dist_factor(2.0);
    set_clip_near_cos(0.3);
    set_clip_boundary_dist_factor(4.0);
    set_clip_boundary_cos(0.3);

    set_consensus_position_dist_factor(1.0);
    set_consensus_normal_dist_factor(3.0);
    set_consensus_jitter_dist_factor(0.01);

    /* These don't really belong under "resolution", but... */
    set_range_data_sigma_factor(4.0);
    set_range_data_min_intensity(0.05);
    set_range_data_horizontal_erode(1);
}


set_zipper_resolution(float res)
{
    ZIPPER_RESOLUTION = res;

    update_edge_length_resolution();
    update_fill_resolution();
    update_align_resolution();
    update_eat_resolution();
    update_clip_resolution();
    update_consensus_resolution();
}


float
get_zipper_resolution()
{
    return ZIPPER_RESOLUTION;
}


