using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SoapOpera {
	class SendProduct {
		string name = "";
		double count = 0;
		double cost = 0;
		string size = "";

		public SendProduct(string Name, double Count, double Cost, string Size) {
			this.name = Name;
			this.count = Count;
			this.cost = Cost;
			this.size = Size;
		}

		public SendProduct(string Name, double Count, double Cost) {
			this.name = Name;
			this.count = Count;
			this.cost = Cost;
			this.size = "";
		}

		public SendProduct() {
			this.name = "";
			this.count = 0;
			this.cost = 0;
			this.size = "";
		}

		public string Name {
			get { return name; }
			set { name = value; }
		}

		public double Count {
			get { return count; }
			set { count = value; }
		}

		public double Cost {
			get { return cost; }
			set { cost = value; }
		}

		public string Size {
			get { return size; }
			set { size = value; }
		}
	}
}
