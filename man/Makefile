TARGETS = Bfield.3.gz \
		  create_angle.3.gz \
		  create_angle_rad.3.gz \
		  create_angle_deg.3.gz \
		  destroy_angle.3.gz \
		  copy_angle.3.gz

all: $(TARGETS)

%.gz: %
	gzip -c $< > $@
