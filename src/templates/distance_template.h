/* special for recalculating delay @# */
 const double diff_x = fabs( x_position[ii] - x_position[target] );
 const double distance_x = fmin( L - diff_x, diff_x );
 const double diff_y = fabs( y_position[ii] - y_position[target] );
 const double distance_y = fmin( L - diff_y, diff_y );

 const double distance = sqrt( distance_x*distance_x + distance_y*distance_y );

 const int delay_in_bins = (int) round( ( synapse_delay + ( distance / conduction_speed[ii] ) ) / dt );
 const int insert_position = ( buffer_position + delay_in_bins ) % buffer_length;

