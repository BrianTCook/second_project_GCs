#!/bin/bash

/usr/bin/expect << EOF
	spawn ./AMUSE_enbid.sh
	expect {
		"Type 0" {
			send -- "1\n"
			}
		}
	}
EOF
