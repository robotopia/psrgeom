TARGETS = Bfield.3.gz \
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
		  set_pulsar.3.gz

all: $(TARGETS)

%.gz: %
	gzip -c $< > $@

clean:
	$(RM) $(TARGETS)
