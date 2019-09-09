using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SoapOpera {
	class SendPeople {
		long id = 0;
		string name = "";
		double allCost = 0;
		List<SendProduct> products;

		public SendPeople(long id, string name, double allCost) {
			this.id = id;
			this.name = name;
			this.allCost = allCost;
			products = new List<SendProduct>();
		}

		public SendPeople() {
			this.id = 0;
			this.name = "";
			this.allCost = 0;
			products = new List<SendProduct>();
		}

		public long Id {
			get { return id; }
			set { id = value; }
		}

		public string Name {
			get { return name; }
			set { name = value; }
		}

		public List<SendProduct> Product {
			get { return products; }
			set { products = value; }
		}

		public double AllCost {
			get { return allCost; }
			set { allCost = value; }
		}
	}
}
