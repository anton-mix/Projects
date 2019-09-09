using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SoapOpera {
	class PhotoInfo {
		double count, price;
		string photoName;

		public PhotoInfo(double count, double price, string photoName) {
			this.count = count;
			this.price = price;
			this.photoName = photoName;
		}

		public double Count {
			get { return count; }
			set { count = value; }
		}

		public double Price {
			get { return price; }
			set { price = value; }
		}

		public string PhotoName {
			get { return photoName; }
			set { photoName = value; }
		}
	}
}
