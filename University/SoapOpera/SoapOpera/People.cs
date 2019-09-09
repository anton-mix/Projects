using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SoapOpera {
	class People {
		string name;
		long id;
		List<List<string>> commentsText;
		List<List<DateTime>> commentsData;
		double count;
		public People(int CountPhoto, string Name, long Id) {
			this.name = Name;
			this.id = Id;
			CommentsText = new List<List<string>>();
			for(int i = 0;i < CountPhoto; i++) {
				CommentsText.Add(new List<string>());
			}
			CommentsData = new List<List<DateTime>>();
			for (int i = 0; i < CountPhoto; i++) {
				CommentsData.Add(new List<DateTime>());
			}
		}

		public string Name {
			get { return name; }
			set { name = value; }
		}

		public long Id {
			get { return id; }
			set { id = value; }
		}

		public List<List<string>> CommentsText {
			get { return commentsText; }
			set { commentsText = value; }
		}

		public List<List<DateTime>> CommentsData {
			get { return commentsData; }
			set { commentsData = value; }
		}
		public double Count {
			get { return count; }
			set { count = value; }
		}
	}
}
