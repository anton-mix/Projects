using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Shapes;
namespace prob1 {
	class MyPoint {
		public  Line LineF;
		public int Z;
		public MyPoint(Line LineF, int Z) {
			this.LineF = LineF;
			this.Z = Z;
		}
	}
}
