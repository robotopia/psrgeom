TARGETS = calc_fields.3.gz \
		  Bstep.3.gz \
		  Bdotrxy.3.gz \
		  footpoint.3.gz \
		  farpoint.3.gz \
		  create_psr_angle.3.gz \
		  create_psr_angle_rad.3.gz \
		  create_psr_angle_deg.3.gz \
		  destroy_psr_angle.3.gz \
		  copy_psr_angle.3.gz \
		  set_psr_angle_rad.3.gz \
		  set_psr_angle_deg.3.gz \
		  set_psr_angle_sin.3.gz \
		  set_psr_angle_cos.3.gz \
		  rotate_about_axis.3.gz \
		  min_phase_diff.3.gz \
		  set_point_xyz.3.gz \
		  set_point_sph.3.gz \
		  set_point_cyl.3.gz \
		  copy_point.3.gz \
		  norm_dot.3.gz \
		  set_pulsar.3.gz \
		  create_pulsar.3.gz \
		  destroy_pulsar.3.gz \
		  set_pulsar_period.3.gz \
		  obs_to_mag_frame.3.gz \
		  mag_to_obs_frame.3.gz \
		  calc_dipole_R.3.gz \
		  dipole_footpoint.3.gz \
		  move_around_cyl.3.gz \
		  print_psrg_header.3.gz \
		  psrgeom.3.gz \
		  psr_cost_lofl.3.gz \
		  psr_cost_los.3.gz \
		  find_emission_point.3.gz \
		  find_approx_emission_point.3.gz \
		  beamangle_to_posangle.3.gz \
		  posangle_to_beamangle.3.gz \
		  line_of_sight.3.gz \
		  psr_fields.7.gz

all: $(TARGETS)

%.gz: %
	gzip -c $< > $@

clean:
	$(RM) $(TARGETS)
