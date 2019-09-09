using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SoapOpera {
	class SendAlbum {
		string name = "";
		List<SendPeople> peoples;
		List<string> photos;
		bool isuseSize = false;
		public SendAlbum(string name) {
			this.name = name;
			peoples = new List<SendPeople>();
			photos = new List<string>();
		}

		public SendAlbum() {
			this.name = "";
			peoples = new List<SendPeople>();
			photos = new List<string>();
		}

		public string Name {
			get { return name; }
			set { name = value; }
		}

		public List<SendPeople> Peoples {
			get { return peoples; }
			set { peoples = value; }
		}

		public List<string> Photos {
			get { return photos; }
			set { photos = value; }
		}

		public bool isUseSize {
			get { return isuseSize; }
			set { isuseSize = value; }
		}
	}
}
