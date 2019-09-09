using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Shapes;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Animation;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;


namespace prob1 {
	class MyStruct {
		//Bild Struct
		Grid MyCanvas;
		int maxNumderFerst = 20,curNumderFerst, rad1 = 15, rad2 = 10, aprLenght,dist = 10;
		double centerXC,centerYC;
		Node[] ferstStruct;
		Random rand = new Random();
		TranslateTransform Trans = new TranslateTransform();
		//Animation Struct
		DoubleAnimation XAnim = new DoubleAnimation(), YAnim = new DoubleAnimation();
		double SpeedPassivAnim = 0.25, SpeedActivAnim = 0.75, offsetXAnim = 0, offsetYAnim = 0,maxPassAnimStat = 50;
		int pastDir = 0,curDir = 0, maxpassAnimRand = 20;
		double transCenterX, transCenterY;
		public MyStruct(Grid MyCanvas,double centerXC, double centerYC) {
			this.MyCanvas = MyCanvas;
			this.centerXC = centerXC;
			this.centerYC = centerYC;
			aprLenght = (int)(((MyCanvas.Height + MyCanvas.Width)/2) * 0.28);
			XAnim.Duration = TimeSpan.FromSeconds(SpeedPassivAnim);
			YAnim.Duration = TimeSpan.FromSeconds(SpeedPassivAnim);
			transCenterX = Trans.X;
			transCenterY = Trans.Y;
		}
		public double TransCenterX {
			get { return transCenterX; }
			set { Trans.X = transCenterX = value; }
		}
		public double TransCenterY {
			get { return transCenterY; }
			set { Trans.Y = transCenterY = value; }
		}
		public void InitFerstStruct() {
			int lineDir;
			double buffX, buffY, dec;
			List<MyPoint> freeVer = new List<MyPoint>();
			freeVer.Add(new MyPoint(null,0));
			curNumderFerst = maxNumderFerst;
			ferstStruct = new Node[curNumderFerst];
			ferstStruct[0] = new Node((int)centerXC,(int)centerYC,rad1,1,MyCanvas,1,Trans);
			for (int i = 0; i < rand.Next(3,5); i++) {
				while (true) {
					if (ferstStruct[0].CountFree == 0)
						break;
					lineDir = ferstStruct[0].NextLine(rand.Next(ferstStruct[0].CountFree));
					if (ferstStruct[0].LengthL[lineDir] == 0) {
						ferstStruct[0].LengthL[(lineDir + 1 == 8) ? 0 : lineDir + 1] = -1;
						ferstStruct[0].LengthL[(lineDir - 1 == -1) ? 7 : lineDir - 1] = -1;
						ferstStruct[0].DrawLine((LineDir)lineDir,rand.Next(aprLenght,(int)(aprLenght * 1.5)),rand.Next(-rad1 * 2,rad1 * 2),rand.Next(-rad1 * 2,rad1 * 2));
						freeVer.Add(new MyPoint(ferstStruct[0].LineF[lineDir],(lineDir + 4 >= 8) ? lineDir - 4 : lineDir + 4));
						break;
					}
				}
			}

			for (int j = 1; j < curNumderFerst; j++) {
				ferstStruct[j] = new Node((int)freeVer[j].LineF.X2,(int)freeVer[j].LineF.Y2,rad1,1,MyCanvas,1,Trans);
				ferstStruct[j].LengthL[freeVer[j].Z] = -2;
				ferstStruct[j].LineF[freeVer[j].Z] = freeVer[j].LineF;
				ferstStruct[j].LengthL[(freeVer[j].Z + 1 == 8) ? 0 : freeVer[j].Z + 1] = -1;
				ferstStruct[j].LengthL[(freeVer[j].Z - 1 == -1) ? 7 : freeVer[j].Z - 1] = -1;
				for (int i = 0; i < rand.Next(1,5); i++) {
					while (true) {
						if (freeVer.Count < curNumderFerst) {
							if (ferstStruct[j].CountFree == 0)
								break;
							lineDir = ferstStruct[j].NextLine(rand.Next(ferstStruct[j].CountFree));
							if (ferstStruct[j].LengthL[lineDir] == 0) {
								ferstStruct[j].DrawLine((LineDir)lineDir,rand.Next(aprLenght,(int)(aprLenght * 1.5)),rand.Next(-rad1 * 2,rad1 * 2),rand.Next(-rad1 * 2,rad1 * 2));
								freeVer.Add(new MyPoint(ferstStruct[j].LineF[lineDir],(lineDir + 4 >= 8) ? lineDir - 4 : lineDir + 4));
								break;
							}
						}
						else
							break;
					}
				}
			}
			for (int i = 0; i < curNumderFerst; i++)
				for (int j = 0; j < curNumderFerst; j++) {
					if (j != i) {
						buffX = Math.Abs(ferstStruct[i].Center.X - ferstStruct[j].Center.X);
						buffY = Math.Abs(ferstStruct[i].Center.Y - ferstStruct[j].Center.Y);
						dec = ferstStruct[i].Rad + ferstStruct[j].Rad + dist;
						if (buffX < dec && buffY < dec) {
							ferstStruct[j].MoveNode((ferstStruct[i].Center.X > ferstStruct[j].Center.X) ? buffX - dec : dec - buffX,(ferstStruct[i].Center.Y > ferstStruct[j].Center.Y) ? dec - buffY : buffY - dec);
						}
					}
				}
		}
		public void deleteStruct() {
			MyCanvas.Children.Clear();
		}

		public void PassivAnimation() {
			while (true) {
				curDir = rand.Next(4);
				if (curDir == pastDir)
					curDir = rand.Next(4);
				else {
					pastDir = curDir;
					break;
				}
			}
			switch (curDir) {
				case 0:
					offsetXAnim = (transCenterX - maxPassAnimStat - rand.Next(maxpassAnimRand)) - Trans.X;
					offsetYAnim = (transCenterY - maxPassAnimStat - rand.Next(maxpassAnimRand)) - Trans.Y;
					break;
				case 1:
					offsetXAnim = (transCenterX + maxPassAnimStat + rand.Next(maxpassAnimRand)) - Trans.X;
					offsetYAnim = (transCenterY - maxPassAnimStat - rand.Next(maxpassAnimRand)) - Trans.Y;
					break;
				case 2:
					offsetXAnim = (transCenterX + maxPassAnimStat + rand.Next(maxpassAnimRand)) - Trans.X;
					offsetYAnim = (transCenterY + maxPassAnimStat + rand.Next(maxpassAnimRand)) - Trans.Y;
					break;
				case 3:
					offsetXAnim = (TransCenterX - maxPassAnimStat - rand.Next(maxpassAnimRand)) - Trans.X;
					offsetYAnim = (TransCenterY + maxPassAnimStat + rand.Next(maxpassAnimRand)) - Trans.Y;
					break;
			}
			XAnim.From = Trans.X;
			YAnim.From = Trans.Y;
			XAnim.To = Trans.X + offsetXAnim;
			YAnim.To = Trans.Y + offsetYAnim;
			Trans.BeginAnimation(TranslateTransform.XProperty,XAnim);
			Trans.BeginAnimation(TranslateTransform.YProperty,YAnim);
		}
		public void ActivAnimation() {

		}

	}
}
