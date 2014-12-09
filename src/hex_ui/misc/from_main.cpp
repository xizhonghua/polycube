int param(int argc, char *argv[])
{
	if(argc < 6) {
		cerr << "param tet zyz wave_length iter_num [More|LBFGS] [theta]" << endl;
		return __LINE__;
	}

	matrix<double> zyz, fixed_frame, aligned, theta;
	matrix<size_t> fixed_frame_idx, aligned_idx;

	tetmesh tm;
	if(tet_mesh_read_from_zjumat(argv[1], &tm.node, &tm.tet))
		return __LINE__;

	cerr << "#node: " << tm.node.size(2)
	   << ", #tet: " << tm.tet.size(2) << endl;
	if(load_from_tet(tm.node, tm.tet,
					 fixed_frame, fixed_frame_idx,
					 aligned, aligned_idx)) {
		cerr << "load fail." << endl;
		return __LINE__;
	}

	ifstream ifs(argv[2], ifstream::binary);
	read_matrix(ifs, zyz);

	cerr << zyz(colon(), colon(0, 6)) << endl;
	cerr << "# of zyz: " << zyz.size(2) << endl;
	matrix<matrix<double> >frame(zyz.size(2));
	for(size_t i = 0; i < frame.size(); ++i) {
		frame[i].resize(3, 3);
		zyz_angle_2_rotation_matrix(zyz(0, i), zyz(1, i), zyz(2, i),
									&frame[i][0]);
		frame[i] = eye<double>(3);
	}

	const double wave_length = wave_length_from_argv(argv[3], tm.node);

	param_opt opt;
	opt.setup_equations(tm.tet, tm.node,
						frame,
						aligned_idx,
						wave_length);

	const int iter_num = atoi(argv[4]);

	if(argc > 6) {
		ifstream ifs(argv[6], ifstream::binary);
		if(!ifs.fail()) {
			read_matrix(ifs, theta);
			if(theta.size(2) != opt.get()->dim_of_x()/3) {
				cerr << "wrong theta file." << endl;
				return __LINE__;
			}
		}
	}
	if(theta.size() == 0) {
		theta = zeros<double>(3, opt.get()->dim_of_x()/3)+3.1415926/4.0;
		auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tm.tet));
		for(size_t fi = 0; fi < fa->faces_.size(); ++fi) {
			matrix<double> ct = zeros<double>(3, 1);
			for(size_t i = 0; i < 3; ++i) {
				ct += tm.node(colon(), fa->faces_[fi][i]);
			}
			ct /= 3;
			// converges to good result below 1e-1, and converges to zero amplitude when over 4e-1
			theta(colon(), fi) = ct*(2*3.1415926)/wave_length;//+rand<double>(3, 1)*2e-1;
//			theta(colon(), fi) = rand<double>(3, 1)*1e1;
		}
		theta = (rand<double>(3, opt.get()->dim_of_x()/3)-0.5);
		cerr << "theta0: " << trans(theta(colon(), 0)) << endl;
	}
	else {
		cout << "get initial value from file." << endl;
		cout << theta(colon(), colon(0, 6));
	}

	matrix<double> residual(opt.get()->dim_of_f());
	opt.solve(theta, residual, iter_num, 1, argv[5]);


	if(argc > 6) {
		ofstream ofs(argv[6], ofstream::binary);
		write_matrix(ofs, theta);
		cout << "theta: " << theta(colon(), colon(0, 6));
		cout << "residual: " << trans(residual(colon(0, 18))) << endl;
//		cout << "residual: " << residual << endl;
	}
	cerr << max(theta) << " " << min(theta) << endl;
	theta = cos(theta);
	cerr << theta(colon(), colon(0, 6)) << endl;
	cerr << max(theta) << " " << min(theta) << endl;
	matrix<double> amp(theta.size(2));
	for(size_t i = 0; i < theta.size(2); ++i)
		amp[i] = theta(0, i)*theta(1, i)*theta(2, i);
	cerr << max(amp) << " " << min(amp) << endl;
	return 0;
}

int param2vtk(int argc, char *argv[])
{
	if(argc < 5) {
		cerr << "param2vtk tet zyz wave_length theta" << endl;
		return __LINE__;
	}

	matrix<double> theta;
	matrix<matrix<double> >frame;

	tetmesh tm;
	if(tet_mesh_read_from_zjumat(argv[1], &tm.node, &tm.tet))
		return __LINE__;

	{
		matrix<double> zyz;
		ifstream ifs(argv[2], ifstream::binary);
		read_matrix(ifs, zyz);
		frame.resize(zyz.size(2));
		for(size_t i = 0; i < frame.size(); ++i) {
			frame[i].resize(3, 3);
			zyz_angle_2_rotation_matrix(zyz(0, i), zyz(1, i), zyz(2, i),
										&frame[i][0]);
			frame[i] = eye<double>(3);
		}
	}

	const double wave_length = wave_length_from_argv(argv[3], tm.node);

	{
		ifstream ifs(argv[4], ifstream::binary);
		cerr << "read theta from " << argv[4] << endl;
		read_matrix(ifs, theta);
	}

	cerr << theta(colon(), colon(0, 6)) << endl;

	auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tm.tet));

	std::vector<float> node_value(tm.node.size(2));
	std::vector<vector<float> > node_values(tm.node.size(2));
	std::vector<size_t> count(tm.node.size(2));
	static const double PI = 3.1415926;

	for(size_t ti = 0; ti < tm.tet.size(2); ++ti) {
		if(ti > 10) continue;
//		cerr << ti << "********************" << endl;
		matrix<double> P = ones<double>(3, 4);
		matrix<size_t> node_idx_(4);
		for(int i = 0; i < 4; ++i) {
			node_idx_[i] = fa->get_face_idx(tm.tet(i, ti), tm.tet((i+1)%4, ti), tm.tet((i+2)%4, ti));
			if(node_idx_[i] >= fa->faces_.size())
				cerr << "wrong face index in fa." << endl;
			matrix<double> face_center = zeros<double>(3, 1);
			for(int j = 0; j < 3; ++j)
				face_center += tm.node(colon(), tm.tet((j+i)%4, ti));
			P(colon(), i) = face_center/3.0;
		}
		static const double PI = 3.1415926;
		matrix<double> edge;
		for(int i = 0; i < 4; ++i) {
			matrix<double> framei = frame[node_idx_[i]];
			for(int j = 0; j < 4; ++j) {
				if(i == j) continue;
				edge = P(colon(), j)-P(colon(), i);
				double amp_ij = 1;
				for(int d = 0; d < 3; ++d) {
					const double dist  = 2*PI*dot(edge, framei(colon(), d))/wave_length;
					amp_ij *= cos(theta(d, node_idx_[i])+dist);
				}
				double amp_jj = 1;
				for(int d = 0; d < 3; ++d) {
					amp_jj *= cos(theta(d, node_idx_[j]));
				}
				const double err = fabs(amp_ij-amp_jj);
//				cerr << "err at theta node: " << amp_ij << " " << amp_jj << " " << err << endl;
			}
			double amp = 1;
			const matrix<double> offset = tm.node(colon(), tm.tet(0, ti))-P(colon(), i);
			for(int d = 0; d < 3; ++d) {
				const double dist = 2*PI*dot(offset, framei(colon(), d))/wave_length;
				amp *= cos(theta(d, node_idx_[i])+dist);
			}
//			cerr << amp << endl;
		}
//		cerr << endl;
	}

	for(size_t fi = 0; fi < fa->faces_.size(); ++fi) { // for face node
		const std::vector<size_t> &facei = fa->faces_[fi];
		matrix<double> theta_pos = zeros<double>(3, 1);
		for(size_t ni = 0; ni < 3; ++ni)
			theta_pos += tm.node(colon(), facei[ni]);
		theta_pos /= 3.0;
		for(size_t ni = 0; ni < 3; ++ni) { // for neighbor tet node
			const matrix<double> offset = tm.node(colon(), facei[ni])-theta_pos;
			double amp = 1;
			for(int d = 0; d < 3; ++d) { // for direction
				const double dist = 2*PI*dot(frame[fi](colon(), d), offset)/wave_length;
				amp *= cos(theta(d, fi)+dist);
			}
			node_value[facei[ni]] += amp;
			node_values[facei[ni]].push_back(amp);
			++count[facei[ni]];
		}
	}
	for(size_t ni = 0; ni < node_value.size(); ++ni) {
		if(node_values[ni].size() != count[ni]) {
			cerr << "strange error." << endl;
		}
		if(ni < 10) {
//			copy(node_values[ni].begin(), node_values[ni].end(), ostream_iterator<float>(cerr, " "));
//			cerr << "\n" << endl;
		}
		node_value[ni] /= count[ni];
	}
//	cerr << *min_element(node_value.begin(), node_value.end()) << " "
//		 << *max_element(node_value.begin(), node_value.end()) << endl;

	cout << "POINT_DATA " << node_value.size() << "\n";
	cout << "SCALARS amplitude float\nLOOKUP_TABLE my_table\n";
	copy(node_value.begin(), node_value.end(), ostream_iterator<float>(cout, "\n"));

	return 0;
}

int param2(int argc, char *argv[])
{
	if(argc < 5) {
		cerr << "param2 tet zyz wave_length iter_num [More|LBFGS] [theta]" << endl;
		return __LINE__;
	}

	matrix<double> zyz, fixed_frame, aligned, theta;
	matrix<size_t> fixed_frame_idx, aligned_idx;

	tetmesh tm;
	if(tet_mesh_read_from_zjumat(argv[1], &tm.node, &tm.tet))
		return __LINE__;

	cerr << "#node: " << tm.node.size(2)
	   << ", #tet: " << tm.tet.size(2) << endl;
	if(load_from_tet(tm.node, tm.tet,
					 fixed_frame, fixed_frame_idx,
					 aligned, aligned_idx)) {
		cerr << "load fail." << endl;
		return __LINE__;
	}

	ifstream ifs(argv[2], ifstream::binary);
	read_matrix(ifs, zyz);

	cerr << zyz(colon(), colon(0, 6)) << endl;
	cerr << "# of zyz: " << zyz.size(2) << endl;
	matrix<matrix<double> >frame(zyz.size(2));
	for(size_t i = 0; i < frame.size(); ++i) {
		frame[i].resize(3, 3);
		zyz_angle_2_rotation_matrix(zyz(0, i), zyz(1, i), zyz(2, i),
									&frame[i][0]);
		frame[i] = eye<double>(3);
	}

	const double wave_length = wave_length_from_argv(argv[3], tm.node);

	param_opt2 opt;
	opt.setup_equations(tm.tet, tm.node,
						frame,
						aligned_idx,
						wave_length);

	const int iter_num = atoi(argv[4]);

	if(argc > 6) {
		ifstream ifs(argv[6], ifstream::binary);
		if(!ifs.fail()) {
			read_matrix(ifs, theta);
			if(theta.size(2) != opt.get()->dim_of_x()/6) {
				cerr << "wrong theta file." << endl;
				return __LINE__;
			}
		}
	}
	if(theta.size() == 0) {
		theta = zeros<double>(6, opt.get()->dim_of_x()/6)+3.1415926/4.0;
		auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tm.tet));
		for(size_t fi = 0; fi < fa->faces_.size(); ++fi) {
			matrix<double> ct = zeros<double>(3, 1);
			for(size_t i = 0; i < 3; ++i) {
				ct += tm.node(colon(), fa->faces_[fi][i]);
			}
			ct /= 3;
			// [theta]: converges to good result below 1e-1, and converges to zero amplitude when over 4e-1
			// [cos(theta), sin(theta)]: good for 1e0, bad for 2e0
			matrix<double> phase = ct*(2*3.1415926)/wave_length;
//			if(fi != 0)
//				phase += rand<double>(3, 1)*2e0;
			for(int d = 0; d < 3; ++d) {
				theta(d*2, fi) = cos(phase[d]);
				theta(d*2+1, fi) = sin(phase[d]);
			}
		}
		theta = (rand<double>(6, opt.get()->dim_of_x()/6)-0.5);
		cerr << "theta0: " << trans(theta(colon(), 0)) << endl;
	}
	else {
		cout << "get initial value from file." << endl;
		cout << theta(colon(), colon(0, 6));
	}

	matrix<double> residual(opt.get()->dim_of_f());
	opt.solve(theta, residual, iter_num, 1, argv[5]);

	if(argc > 6) {
		ofstream ofs(argv[6], ofstream::binary);
		write_matrix(ofs, theta);
		cout << "theta: " << theta(colon(), colon(0, 6));
		cout << "residual: " << trans(residual(colon(0, 18))) << endl;
		matrix<double> fabs_residual = fabs(residual);
		cout << "max residual is: " << max(fabs_residual)
			 << " at: " << max_element(fabs_residual.begin(), fabs_residual.end())-fabs_residual.begin()
			 << endl;

		cout << "fix residual: " << residual(colon(residual.size()-6, residual.size()-1)) << endl;
	}
	cerr << theta(colon(), colon(0, 6)) << endl;
	cerr << max(theta) << " " << min(theta) << endl;
	matrix<double> amp(theta.size(2));
	for(size_t i = 0; i < theta.size(2); ++i)
		amp[i] = theta(0, i)*theta(2, i)*theta(4, i);
	cerr << max(amp) << " " << min(amp) << endl;
	return 0;
}

int param2vtk2(int argc, char *argv[])
{
	if(argc < 5) {
		cerr << "param2vtk2 tet zyz wave_length theta" << endl;
		return __LINE__;
	}

	matrix<double> theta;
	matrix<matrix<double> >frame;

	tetmesh tm;
	if(tet_mesh_read_from_zjumat(argv[1], &tm.node, &tm.tet))
		return __LINE__;

	{
		matrix<double> zyz;
		ifstream ifs(argv[2], ifstream::binary);
		read_matrix(ifs, zyz);
		frame.resize(zyz.size(2));
		for(size_t i = 0; i < frame.size(); ++i) {
			frame[i].resize(3, 3);
			zyz_angle_2_rotation_matrix(zyz(0, i), zyz(1, i), zyz(2, i),
										&frame[i][0]);
			frame[i] = eye<double>(3);
		}
	}

	const double wave_length = wave_length_from_argv(argv[3], tm.node);

	{
		ifstream ifs(argv[4], ifstream::binary);
		read_matrix(ifs, theta);
	}

//	cerr << theta(colon(), colon(0, 6)) << endl;

	auto_ptr<face2tet_adjacent> fa(face2tet_adjacent::create(tm.tet));

	std::vector<float> node_value(tm.node.size(2));
	std::vector<vector<float> > node_values(tm.node.size(2));
	std::vector<size_t> count(tm.node.size(2));
	static const double PI = 3.1415926;

	for(size_t ti = 0; ti < tm.tet.size(2); ++ti) {
		if(ti > 10) continue;
//		cerr << ti << "********************" << endl;
		matrix<double> P = ones<double>(3, 4);
		matrix<size_t> node_idx_(4);
		for(int i = 0; i < 4; ++i) {
			node_idx_[i] = fa->get_face_idx(tm.tet(i, ti), tm.tet((i+1)%4, ti), tm.tet((i+2)%4, ti));
			if(node_idx_[i] >= fa->faces_.size())
				cerr << "wrong face index in fa." << endl;
			matrix<double> face_center = zeros<double>(3, 1);
			for(int j = 0; j < 3; ++j)
				face_center += tm.node(colon(), tm.tet((j+i)%4, ti));
			P(colon(), i) = face_center/3.0;
		}
		static const double PI = 3.1415926;
		matrix<double> edge;
		for(int i = 0; i < 4; ++i) {
			matrix<double> framei = frame[node_idx_[i]];
			for(int j = 0; j < 4; ++j) {
				if(i == j) continue;
				edge = P(colon(), j)-P(colon(), i);
				double amp_ij = 1;
				for(int d = 0; d < 3; ++d) {
					const double dist  = 2*PI*dot(edge, framei(colon(), d))/wave_length;
					amp_ij *= theta(d*2, node_idx_[i])*cos(dist)
						-theta(d*2+1, node_idx_[i])*sin(dist);
				}
				double amp_jj = 1;
				for(int d = 0; d < 3; ++d) {
					amp_jj *= theta(d*2, node_idx_[j]);
				}
				const double err = fabs(amp_ij-amp_jj);
//				cerr << "err at theta node: " << amp_ij << " " << amp_jj << " " << err << endl;
			}
			double amp = 1;
			const matrix<double> offset = tm.node(colon(), tm.tet(0, ti))-P(colon(), i);
			for(int d = 0; d < 3; ++d) {
				const double dist = 2*PI*dot(offset, framei(colon(), d))/wave_length;
				amp *= theta(d*2, node_idx_[i])*cos(dist)
					-theta(d*2+1, node_idx_[i])*sin(dist);
			}
//			cerr << amp << endl;
		}
//		cerr << endl;
	}

	for(size_t fi = 0; fi < fa->faces_.size(); ++fi) { // for face node
		const std::vector<size_t> &facei = fa->faces_[fi];
		matrix<double> theta_pos = zeros<double>(3, 1);
		for(size_t ni = 0; ni < 3; ++ni)
			theta_pos += tm.node(colon(), facei[ni]);
		theta_pos /= 3.0;
		for(size_t ni = 0; ni < 3; ++ni) { // for neighbor tet node
			const matrix<double> offset = tm.node(colon(), facei[ni])-theta_pos;
			double amp = 1;
			for(int d = 0; d < 3; ++d) { // for direction
				const double dist = 2*PI*dot(frame[fi](colon(), d), offset)/wave_length;
				amp *= theta(d*2, fi)*cos(dist)
					-theta(d*2+1, fi)*sin(dist);
			}
			node_value[facei[ni]] += amp;
			node_values[facei[ni]].push_back(amp);
			++count[facei[ni]];
		}
	}
	for(size_t ni = 0; ni < node_value.size(); ++ni) {
		if(node_values[ni].size() != count[ni]) {
			cerr << "strange error." << endl;
		}
		if(ni < 10) {
//			copy(node_values[ni].begin(), node_values[ni].end(), ostream_iterator<float>(cerr, " "));
//			cerr << "\n" << endl;
		}
		node_value[ni] /= count[ni];
	}
//	cerr << *min_element(node_value.begin(), node_value.end()) << " "
//		 << *max_element(node_value.begin(), node_value.end()) << endl;

	cout << "POINT_DATA " << node_value.size() << "\n";
	cout << "SCALARS amplitude float\nLOOKUP_TABLE my_table\n";
	copy(node_value.begin(), node_value.end(), ostream_iterator<float>(cout, "\n"));

	return 0;
}
