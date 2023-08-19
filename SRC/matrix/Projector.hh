Project::Deviator
Project::Spherical

// from DruckerPrager.cpp

    // 2nd order Identity Tensor
    mI1.Zero();
    mI1(0) = 1;
    mI1(1) = 1;
    mI1(2) = 1;

    // 4th order Volumetric Tensor
    // IIvol = I1 tensor I1
    mIIvol.Zero();
    mIIvol(0,0) = 1;
    mIIvol(0,1) = 1;
    mIIvol(0,2) = 1;
    mIIvol(1,0) = 1;
    mIIvol(1,1) = 1;
    mIIvol(1,2) = 1;
    mIIvol(2,0) = 1;
    mIIvol(2,1) = 1;
    mIIvol(2,2) = 1;

    // 4th order Deviatoric Tensor
    // Note:  this is the contravariant form!
    //        usable for s^a = 2G * IIdev^ab * epsilon_b
    // (Need a different form for s^a = IIdev ^a_b * sigma^a)
    mIIdev.Zero();
    mIIdev(0,0) = two3;
    mIIdev(0,1) = -one3;
    mIIdev(0,2) = -one3;
    mIIdev(1,0) = -one3;
    mIIdev(1,1) = two3;
    mIIdev(1,2) = -one3;
    mIIdev(2,0) = -one3;
    mIIdev(2,1) = -one3;
    mIIdev(2,2) = two3;
    mIIdev(3,3) = 0.5;
    mIIdev(4,4) = 0.5;
    mIIdev(5,5) = 0.5;

// From ManzariDafalias
Vector              ManzariDafalias::mI1(6);
Matrix              ManzariDafalias::mIIco(6,6);
Matrix              ManzariDafalias::mIIcon(6,6);
Matrix              ManzariDafalias::mIImix(6,6);
Matrix              ManzariDafalias::mIIvol(6,6);
Matrix              ManzariDafalias::mIIdevCon(6,6);
Matrix              ManzariDafalias::mIIdevMix(6,6);
Matrix              ManzariDafalias::mIIdevCo(6,6);

		initTensors() {
			// 2nd order identity tensor
			mI1.Zero();
			mI1(0) = 1.0;
			mI1(1) = 1.0;
			mI1(2) = 1.0;
			// 4th order mixed variant identity tensor
			mIImix.Zero();
			for (int i = 0; i<6; i++) {
				mIImix(i,i) = 1.0;
			}
			// 4th order covariant identity tensor
			mIIco = mIImix;
			mIIco(3,3) = 2.0;
			mIIco(4,4) = 2.0;
			mIIco(5,5) = 2.0;
			// 4th order contravariant identity tensor
			mIIcon = mIImix;
			mIIcon(3,3) = 0.5;
			mIIcon(4,4) = 0.5;
			mIIcon(5,5) = 0.5;
			// 4th order volumetric tensor, IIvol = I1 tensor I1
			mIIvol.Zero();
			for (int i = 0; i<3; i++) {
				mIIvol(i,0) = 1.0;
				mIIvol(i,1) = 1.0;
				mIIvol(i,2) = 1.0;
			}
			// 4th order contravariant deviatoric tensor
			mIIdevCon = mIIcon - one3*mIIvol;
			// 4th order covariant deviatoric tensor
			mIIdevCo  = mIIco  - one3*mIIvol;
			// 4th order mixed variant deviatoric tensor
			mIIdevMix = mIImix - one3*mIIvol;
		}
	} initTensorOps;
